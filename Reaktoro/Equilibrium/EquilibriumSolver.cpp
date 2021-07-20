// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
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
namespace {

/// Return an EquilibriumSpecs object that represent the specifications of classic Gibbs energy minimization problem.
auto defaultEquilibriumSpecs(const ChemicalSystem& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    return specs;
}

} // namespace

struct EquilibriumSolver::Impl
{
    /// The chemical system associated with this equilibrium solver.
    const ChemicalSystem system;

    /// The chemical equilibrium specifications associated with this equilibrium solver.
    const EquilibriumSpecs specs;

    /// The dimensions of the variables and constraints in the equilibrium specifications.
    const EquilibriumDims dims;

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

    /// Construct a Impl instance with given EquilibriumConditions object.
    Impl(const EquilibriumSpecs& specs)
    : system(specs.system()), specs(specs), dims(specs), setup(specs)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    /// Set the options of the equilibrium solver.
    auto setOptions(const EquilibriumOptions& opts) -> void
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
            // Use `n` instead of `x` to name the variables
            options.optima.output.xprefix = "n";

            // Define some auxiliary references to the variables names
            auto& xnames = options.optima.output.xnames;

            // Initialize the names of the variables corresponding to the species
            for(auto species : system.species())
                xnames.push_back(species.name());

            // Initialize the names of the variables corresponding to the implicit titrants
            for(auto titrant : specs.namesTitrantsImplicit())
                xnames.push_back(titrant);
        }

        // Pass along the options used for the calculation to Optima::Solver object
        optsolver.setOptions(options.optima);
    }

    /// Update the initial chemical state before a new equilibrium calculation.
    auto updateInitialState(ChemicalState& state0, const EquilibriumConditions& conditions)
    {
        const auto n0 = conditions.initialSpeciesAmounts();
        const auto uninitialized = state0.speciesAmounts().maxCoeff() == 0.0;
        const auto has_n0 = n0.size() > 0;

        errorif(uninitialized && !has_n0, "Provide an initialized chemical state as initial guess for the calculation when initial component amounts are given.");

        // Change species amounts only if all amounts are zero
        if(uninitialized && has_n0)
            state0.setSpeciesAmounts(n0);
    }

    /// Update the optimization problem before a new equilibrium calculation.
    auto updateOptProblem(const ChemicalState& state0, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions)
    {
        // The input variables for the equilibrium calculation
        const VectorXr w = conditions.inputValues();

        // Create the Optima::Dims object with dimension info of the optimization problem
        optdims = Optima::Dims();
        optdims.x  = dims.Nx;
        optdims.p  = dims.Np;
        optdims.be = dims.Nb;
        optdims.c  = dims.Nw + dims.Nb; // c = (w, b) where w are the input variables and b are the amounts of components

        // Recreate a new Optima::Problem problem (TODO: Avoid recreation of Optima::Problem object for each equilibrium calculation)
        optproblem = Optima::Problem(optdims);

        // Set the objective function in the Optima::Problem object
        optproblem.f = [=](Optima::ObjectiveResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ObjectiveOptions opts) mutable
        {
            res.f = setup.evalObjectiveValue(x, p, w);
            res.fx = setup.evalObjectiveGradX(x, p, w);

            if(opts.eval.fxc) // if true, this means we are now at the step of computing sensitivity derivatives
                setup.assembleChemicalPropsJacobianBegin(); // start recording the derivatives of the chemical properties with respect to n, p, w

            if(opts.eval.fxx)
                res.fxx = setup.evalObjectiveHessianX(x, p, w); // TODO: Implement diagonal approximation mode for Hessian matrix in EquilibriumSolver.

            if(opts.eval.fxp)
                res.fxp = setup.evalObjectiveHessianP(x, p, w);

            if(opts.eval.fxc)
                res.fxc = setup.evalObjectiveHessianW(x, p, w);

            if(opts.eval.fxc)
                setup.assembleChemicalPropsJacobianEnd(); // finalize recording of the derivatives of the chemical properties with respect to n, p, w

            res.succeeded = true;
        };

        // Set the external constraint function in the Optima::Problem object
        optproblem.v = [=](Optima::ConstraintResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ConstraintOptions opts) mutable
        {
            res.val = setup.evalEquationConstraints(x, p, w);

            if(opts.eval.ddx)
                res.ddx = setup.evalEquationConstraintsGradX(x, p, w);

            if(opts.eval.ddp)
                res.ddp = setup.evalEquationConstraintsGradP(x, p, w);

            if(opts.eval.ddc)
                res.ddc = setup.evalEquationConstraintsGradW(x, p, w);

            res.succeeded = true;
        };

        // Set the coefficient matrices Aex and Aep of the linear equality constraints
        optproblem.Aex = setup.assembleMatrixAex();
        optproblem.Aep = setup.assembleMatrixAep();

        /// Set the right-hand side vector be of the linear equality constraints.
        optproblem.be = setup.assembleVectorBe(conditions);

        // Set the lower bounds of the species amounts
        optproblem.xlower = setup.assembleLowerBoundsVector(restrictions, state0);
        optproblem.xupper = setup.assembleUpperBoundsVector(restrictions, state0);

        // Set the values of the input variables for sensitivity derivatives (due to the use of Param, a wrapper to a shared pointer, the actual values of c here are not important, because the Param objects are embedded in the models)
        optproblem.c = zeros(optdims.c);

        // Set the Jacobian matrix d(be)/dc = [d(be)/dw d(be)/db]
        // The left Nw x Nb block is zero. The right Nb x Nb block is identity!
        optproblem.bec.setZero();
        optproblem.bec.rightCols(dims.Nb).diagonal().setOnes();
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateOptState(const ChemicalState& state0)
    {
        optstate = state0.equilibrium().optimaState();

        if(optstate.dims.x != dims.Nx)
            optstate = Optima::State(optdims);

        optstate.x.head(dims.Nn) = state0.speciesAmounts();

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
    auto updateChemicalState(ChemicalState& state, const EquilibriumConditions& conditions)
    {
        const auto T = setup.chemicalProps().temperature();
        const auto P = setup.chemicalProps().pressure();
        state.setTemperature(T);
        state.setPressure(P);
        state.setSpeciesAmounts(optstate.x.head(dims.Nn));
        state.equilibrium().setInputNames(conditions.inputNames());
        state.equilibrium().setInputValues(conditions.inputValues());
        state.equilibrium().setInitialComponentAmounts(optproblem.be);
        state.equilibrium().setOptimaState(optstate);
    }

    /// Update the equilibrium sensitivity object with computed optimization sensitivity.
    auto updateEquilibriumSensitivity(EquilibriumSensitivity& sensitivity)
    {
        const auto& Nn = dims.Nn;
        const auto& Nb = dims.Nb;
        const auto& Nq = dims.Nq;
        const auto& Nw = dims.Nw;
        const auto& xc = optsensitivity.xc;
        const auto& pc = optsensitivity.pc;

        const auto& props = setup.equilibriumProps();

        const auto dndw = xc.topLeftCorner(Nn, Nw);
        const auto dqdw = xc.bottomLeftCorner(Nq, Nw);
        const auto dpdw = pc.leftCols(Nw);
        const auto dndb = xc.topRightCorner(Nn, Nb);
        const auto dqdb = xc.bottomRightCorner(Nq, Nb);
        const auto dpdb = pc.rightCols(Nb);
        const auto dudn = props.dudn();
        const auto dudp = props.dudp();
        const auto dudw = props.dudw();

        sensitivity.initialize(specs);
        sensitivity.dndw(dndw);
        sensitivity.dqdw(dqdw);
        sensitivity.dpdw(dpdw);
        sensitivity.dndb(dndb);
        sensitivity.dqdb(dqdb);
        sensitivity.dpdb(dpdb);
        sensitivity.dudw(dudw + dudn*dndw + dudp*dpdw);
        sensitivity.dudb(dudn*dndb + dudp*dpdb);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state0, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
    {
        EquilibriumConditions conditions(specs);
        conditions.temperature(state0.temperature());
        conditions.pressure(state0.pressure());
        conditions.startWithState(state0);
        return solve(state0, conditions, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, const EquilibriumConditions& conditions) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state0, conditions, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
    {
        EquilibriumResult eqresult;

        updateInitialState(state0, conditions);
        updateOptProblem(state0, conditions, restrictions);
        updateOptState(state0);

        eqresult.optima = optsolver.solve(optproblem, optstate);

        updateChemicalState(state0, conditions);

        return eqresult;
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, EquilibriumSensitivity& sensitivity) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state0, sensitivity, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, EquilibriumSensitivity& sensitivity, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
    {
        EquilibriumConditions conditions(specs);
        conditions.temperature(state0.temperature());
        conditions.pressure(state0.pressure());
        conditions.startWithState(state0);
        return solve(state0, sensitivity, conditions, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions) -> EquilibriumResult
    {
        EquilibriumRestrictions restrictions(system);
        return solve(state0, sensitivity, conditions, restrictions);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
    {
        EquilibriumResult eqresult;

        updateInitialState(state0, conditions);
        updateOptProblem(state0, conditions, restrictions);
        updateOptState(state0);

        eqresult.optima = optsolver.solve(optproblem, optstate, optsensitivity);

        updateChemicalState(state0, conditions);
        updateEquilibriumSensitivity(sensitivity);

        return eqresult;
    }
};

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(defaultEquilibriumSpecs(system)))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumSpecs& specs)
: pimpl(new Impl(specs))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSolver::~EquilibriumSolver()
{}

auto EquilibriumSolver::operator=(EquilibriumSolver other) -> EquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->solve(state);
}

auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumConditions& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions);
}

auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, conditions, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, restrictions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions);
}

auto EquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions) -> EquilibriumResult
{
    return pimpl->solve(state, sensitivity, conditions, restrictions);
}

} // namespace Reaktoro
