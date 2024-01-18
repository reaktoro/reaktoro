// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "EquilibriumSolver.hpp"

// Optima includes
#include <Optima/Options.hpp>
#include <Optima/Problem.hpp>
#include <Optima/Result.hpp>
#include <Optima/Sensitivity.hpp>
#include <Optima/Solver.hpp>
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/ArrayStream.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Warnings.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSetup.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

auto const EQUILIBRIUM_FAILURE_MESSAGE = R"(The chemical equilibrium/kinetics calculation did not converge.
This can occur due to various factors, such as:
  1. Infeasible Problem Formulation: The problem you have defined lacks a feasible solution within the given constraints and assumptions.
  2. Out-of-Range Thermodynamic Conditions: The conditions you have specified for the system may fall outside the valid range supported by the thermodynamic models used.
  3. Numerical Instabilities: Convergence issues may arise from numerical problems during the execution of the chemical/kinetic equilibrium algorithm. Consider reporting the issue with a minimal reproducible example if you believe the algorithm is responsible for this issue.
Disable this warning message with Warnings.disable(906) in Python and Warnings::disable(906) in C++.)";

struct EquilibriumSolver::Impl
{
    /// The chemical system associated with this equilibrium solver.
    const ChemicalSystem system;

    /// The chemical equilibrium specifications associated with this equilibrium solver.
    const EquilibriumSpecs specs;

    /// The dimensions of the variables and constraints in the equilibrium specifications.
    const EquilibriumDims dims;

    /// The auxiliary equilibrium conditions used whenever none are given in the solve methods.
    const EquilibriumConditions xconditions;

    /// The auxiliary equilibrium restrictions used whenever none are given in the solve methods.
    const EquilibriumRestrictions xrestrictions;

    /// The auxiliary amounts of conservative components used whenever none are given in the solve methods.
    const ArrayXd xc0;

    // The equilibrium problem setup for the equilibrium solver.
    EquilibriumSetup setup;

    /// The options of the equilibrium solver.
    EquilibriumOptions options;

    /// The dimensions of the variables and constraints in the optimization problem.
    Optima::Dims optdims;

    /// The optimization problem to be configured for a chemical equilibrium calculation.
    Optima::Problem optproblem;

    /// The optimization state of the calculation.
    Optima::State optstate;

    /// The optimization sensitivity of the calculation.
    Optima::Sensitivity optsensitivity;

    /// The solver for the optimization calculations.
    Optima::Solver optsolver;

    /// The result of the equilibrium calculation
    EquilibriumResult result;

    /// The array stream used to clean up autodiff seed values from the last ChemicalProps update step.
    ArrayStream<double> stream;

    /// Construct a Impl instance with given EquilibriumConditions object.
    Impl(EquilibriumSpecs const& specs)
    : system(specs.system()), specs(specs), dims(specs), xconditions(specs), xrestrictions(system), setup(specs)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    /// Set the options of the equilibrium solver.
    auto setOptions(EquilibriumOptions const& opts) -> void
    {
        // Update the options of the equilibrium calculation
        options = opts;

        // Pass along the options used for the calculation to EquilibriumSetup object
        setup.setOptions(options);

        // Ensure some options have proper values
        error(options.epsilon <= 0, "EquilibriumOptions::epsilon cannot be zero or negative.");

        // Initialize the names of the primal and dual variables
        if(options.optima.output.active)
        {
            // Define some auxiliary references to the variables names
            auto& xnames = options.optima.output.xnames;
            auto& ynames = options.optima.output.ynames;
            auto& pnames = options.optima.output.pnames;

            // Clear the names of the variables
            xnames.clear();
            ynames.clear();
            pnames.clear();

            // Initialize the names of the x variables corresponding to species amounts
            for(auto species : system.species())
                xnames.push_back(species.name());

            // Initialize the names of the x variables corresponding to implicit titrant amounts
            for(auto qname : specs.namesControlVariablesQ())
                xnames.push_back(qname);

            // Initialize the names of the y Lagrange multipliers corresponding to the amounts of conservative components
            for(auto component : specs.namesConservativeComponents())
                ynames.push_back(component);

            // Initialize the names of the p variables corresponding to the explicit titrants, temperature and/or pressure, other added unknowns
            for(auto pname : specs.namesControlVariablesP())
                pnames.push_back(pname);
        }

        // Pass along the options used for the calculation to Optima::Solver object
        optsolver.setOptions(options.optima);
    }

    /// Update the optimization problem before a new equilibrium calculation.
    auto updateOptProblem(ChemicalState const& state0, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions)
    {
        // The input variables for the equilibrium calculation
        const VectorXr w = conditions.inputValuesGetOrCompute(state0);

        // Create the Optima::Dims object with dimension info of the optimization problem
        optdims = Optima::Dims();
        optdims.x  = dims.Nx;
        optdims.p  = dims.Np;
        optdims.be = dims.Nc;
        optdims.c  = dims.Nw + dims.Nc; // c' = (w, c) where w are the input variables and c are the amounts of components

        // Recreate a new Optima::Problem problem (TODO: Avoid recreation of Optima::Problem object for each equilibrium calculation)
        optproblem = Optima::Problem(optdims);

        // Set the resources function in the Optima::Problem object
        optproblem.r = [=](VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ObjectiveOptions fopts, Optima::ConstraintOptions hopts, Optima::ConstraintOptions vopts) mutable
        {
            setup.update(x, p, w);

            if(fopts.eval.fxc || vopts.eval.ddc)
                setup.assembleChemicalPropsJacobianBegin();

            if(fopts.eval.fxx || vopts.eval.ddx)
                setup.updateGradX(fopts.ibasicvars);
            if(fopts.eval.fxp || vopts.eval.ddp)
                setup.updateGradP();
            if(fopts.eval.fxc || vopts.eval.ddc)
                setup.updateGradW();

            if(fopts.eval.fxc || vopts.eval.ddc)
                setup.assembleChemicalPropsJacobianEnd();
        };

        // Set the objective function in the Optima::Problem object
        optproblem.f = [=](Optima::ObjectiveResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ObjectiveOptions opts) mutable
        {
            res.f = setup.getGibbsEnergy();
            res.fx = setup.getGibbsGradX();

            if(opts.eval.fxx)
                res.fxx = setup.getGibbsHessianX();

            if(opts.eval.fxp)
                res.fxp = setup.getGibbsHessianP();

            if(opts.eval.fxc)
                res.fxc = setup.getGibbsHessianC();

            res.fxx4basicvars = setup.usingPartiallyExactDerivatives();
            res.diagfxx = setup.usingDiagonalApproxDerivatives();

            res.succeeded = true;
        };

        // Set the external constraint function in the Optima::Problem object
        optproblem.v = [=](Optima::ConstraintResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ConstraintOptions opts) mutable
        {
            res.val = setup.getConstraintResiduals();

            if(opts.eval.ddx)
                res.ddx = setup.getConstraintResidualsGradX();

            if(opts.eval.ddp)
                res.ddp = setup.getConstraintResidualsGradP();

            if(opts.eval.ddc)
                res.ddc = setup.getConstraintResidualsGradC();

            res.succeeded = true;
        };

        // Set the coefficient matrices Aex and Aep of the linear equality constraints
        optproblem.Aex = setup.Aex();
        optproblem.Aep = setup.Aep();

        /// Set the right-hand side vector be of the linear equality constraints.
        optproblem.be = conditions.initialComponentAmountsGetOrCompute(state0);

        // Set the lower and upper bounds of the species amounts
        optproblem.xlower = setup.assembleLowerBoundsVector(restrictions, state0);
        optproblem.xupper = setup.assembleUpperBoundsVector(restrictions, state0);

        // Set the lower and upper bounds of the *p* control variables
        optproblem.plower = conditions.lowerBoundsControlVariablesP();
        optproblem.pupper = conditions.upperBoundsControlVariablesP();

        // Set the values of the input variables for sensitivity derivatives
        optproblem.c = zeros(optdims.c);

        // Set the Jacobian matrix d(be)/dc = [d(be)/dw d(be)/db]
        // The left Nw x Nb block is zero. The right Nb x Nb block is identity!
        optproblem.bec.setZero();
        optproblem.bec.rightCols(dims.Nc).diagonal().setOnes();
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateOptState(ChemicalState const& state0)
    {
        // Initialize optstate with that stored in state0 (note state0 may have empty Optima::State object!)
        optstate = state0.equilibrium().optimaState();

        // In case optstate corresponds to an equilibrium problem of different structure, initialize it with a clean slate
        if(optstate.dims.x != dims.Nx || optstate.dims.p != dims.Np || optstate.dims.be != dims.Nc || optstate.dims.c != dims.Nw + dims.Nc)  // TODO: Replace this by a code that represents the EquilibriumSpecs object used for the previous calculation. Consider a dictionary of saved optstates and corresponding EquilibriumSpecs objects in case the same ChemicalState object is used within different solvers.
            optstate = Optima::State(optdims);

        // Overwrite n in x = (n, q) with species amounts from the chemical state
        optstate.x.head(dims.Nn) = state0.speciesAmounts();

        // Overwrite delta variables in q with zeros (i.e., the amount of an implicit titrant to add/remove)
        optstate.x.tail(dims.Nq).fill(0.0);

        // Overwrite delta variables in p with zeros (i.e., the amount of an explicit titrant to add/remove, temperature and/or pressure increase/decrease)
        optstate.p.fill(0.0);

        // TODO: Instead of using T and P as unknown, use dT and dP, so this is block of code is not needed!
        if(specs.isTemperatureUnknown())
        {
            optstate.p[0] = state0.temperature();
            if(specs.isPressureUnknown())
                optstate.p[1] = state0.pressure();
        }
        else if(specs.isPressureUnknown())
            optstate.p[0] = state0.pressure();
    }

    /// Update the chemical state object with computed optimization state.
    auto updateChemicalState(ChemicalState& state, EquilibriumConditions const& conditions)
    {
        // Update the ChemicalProps object in state
        auto& props = state.props();
        props = setup.chemicalProps();

        // TODO: In Optima, make sure check for convergence does not compute
        // any derivatives. Use F.updateSkipJacobian(u) instead of F.update(u)
        // in method MasterSolver::Impl::stepping. Once this is implemented,
        // there will be no need for this method, because the chemical
        // properties will be clean of derivatives (i.e., autodiff seed values
        // will be zero). Once this is done, the next step of cleaning up such
        // seed values can be removed.

        // Make sure the derivative information in the underlying chemical
        // properties of the system are zeroed out!
        ArrayStream<double> stream;
        props.serialize(stream);
        props.deserialize(stream);

        // Update other state variables in the ChemicalState object
        state.setTemperature(props.temperature());
        state.setPressure(props.pressure());
        state.setSpeciesAmounts(optstate.x.head(dims.Nn));
        state.equilibrium().setNamesInputVariables(specs.namesInputs());
        state.equilibrium().setNamesControlVariablesP(specs.namesControlVariablesP());
        state.equilibrium().setNamesControlVariablesQ(specs.namesControlVariablesQ());
        state.equilibrium().setInputVariables(conditions.inputValues());
        state.equilibrium().setInitialComponentAmounts(optproblem.be);
        state.equilibrium().setOptimaState(optstate);
    }

    /// Update the equilibrium sensitivity object with computed optimization sensitivity.
    auto updateEquilibriumSensitivity(EquilibriumSensitivity& sensitivity)
    {
        auto const& Nn = dims.Nn;
        auto const& Nc = dims.Nc;
        auto const& Nq = dims.Nq;
        auto const& Nw = dims.Nw;
        auto const& xc = optsensitivity.xc;
        auto const& pc = optsensitivity.pc;

        auto const& props = setup.equilibriumProps();

        const auto dndw = xc.topLeftCorner(Nn, Nw);
        const auto dqdw = xc.bottomLeftCorner(Nq, Nw);
        const auto dpdw = pc.leftCols(Nw);
        const auto dndc = xc.topRightCorner(Nn, Nc);
        const auto dqdc = xc.bottomRightCorner(Nq, Nc);
        const auto dpdc = pc.rightCols(Nc);
        const auto dudn = props.dudn();
        const auto dudp = props.dudp();
        const auto dudw = props.dudw();

        sensitivity.initialize(specs);
        sensitivity.dndw(dndw);
        sensitivity.dqdw(dqdw);
        sensitivity.dpdw(dpdw);
        sensitivity.dndc(dndc);
        sensitivity.dqdc(dqdc);
        sensitivity.dpdc(dpdc);
        sensitivity.dudw(dudw + dudn*dndw + dudp*dpdw);
        sensitivity.dudc(dudn*dndc + dudp*dpdc);
    }

    auto solve(ChemicalState& state) -> EquilibriumResult
    {
        return solve(state, xconditions, xrestrictions);
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return solve(state, xconditions, restrictions);
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        return solve(state, conditions, xrestrictions);
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        updateOptProblem(state, conditions, restrictions);
        updateOptState(state);

        const auto optstatebkp = optstate;

        result.optima = optsolver.solve(optproblem, optstate);

        if(!result.optima.succeeded)
        {
            auto optionsbkp = options;
            options.optima.backtracksearch.apply_min_max_fix_and_accept = !options.optima.backtracksearch.apply_min_max_fix_and_accept;
            setOptions(options);
            optstate = optstatebkp;
            result.optima = optsolver.solve(optproblem, optstate);
            options = optionsbkp;
            setOptions(options);
        }

        warningif(!result.optima.succeeded && Warnings::isEnabled(906), EQUILIBRIUM_FAILURE_MESSAGE);

        updateChemicalState(state, conditions);

        return result;
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> EquilibriumResult
    {
        return solve(state, sensitivity, xconditions, xrestrictions);
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        return solve(state, sensitivity, xconditions, restrictions);
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> EquilibriumResult
    {
        return solve(state, sensitivity, conditions, xrestrictions);
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
    {
        EquilibriumResult result;

        updateOptProblem(state, conditions, restrictions);
        updateOptState(state);

        result.optima = optsolver.solve(optproblem, optstate, optsensitivity);

        updateChemicalState(state, conditions);
        updateEquilibriumSensitivity(sensitivity);

        return result;
    }
};

EquilibriumSolver::EquilibriumSolver(ChemicalSystem const& system)
: pimpl(new Impl(EquilibriumSpecs::TP(system)))
{}

EquilibriumSolver::EquilibriumSolver(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

EquilibriumSolver::EquilibriumSolver(EquilibriumSolver const& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSolver::~EquilibriumSolver()
{}

auto EquilibriumSolver::operator=(EquilibriumSolver other) -> EquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->solve(state);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions, restrictions);
}

auto EquilibriumSolver::setOptions(EquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
