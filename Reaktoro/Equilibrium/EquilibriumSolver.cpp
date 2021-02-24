// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

    /// The auxiliary vector to store the amounts of the species.
    ArrayXd n0;

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

    /// Update the optimization problem before a new equilibrium calculation.
    auto updateOptProblem(ChemicalState& state0, const EquilibriumConditions& conditions, const EquilibriumRestrictions& restrictions)
    {
        // Auxiliary data
        const auto params = conditions.params();

        // Create the Optima::Dims object with dimension info of the optimization problem
        optdims = Optima::Dims();
        optdims.x  = dims.Nx;
        optdims.p  = dims.Np;
        optdims.be = dims.Nb;

        // Recreate a new Optima::Problem problem (TODO: Avoid recreation of Optima::Problem object for each equilibrium calculation)
        optproblem = Optima::Problem(optdims);

        // Set the objective function in the Optima::Problem object
        optproblem.f = [=](Optima::ObjectiveResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ObjectiveOptions opts) mutable
        {
            res.f = setup.evalObjectiveValue(x, p, params);
            res.fx = setup.evalObjectiveGradX(x, p, params);

            if(opts.eval.fxx)
                res.fxx = setup.evalObjectiveHessianX(x, p, params); // TODO: Implement diagonal approximation mode for Hessian matrix in EquilibriumSolver.

            if(opts.eval.fxp)
                res.fxp = setup.evalObjectiveHessianP(x, p, params);

            if(opts.eval.fxc)
                res.fxc = setup.evalObjectiveHessianParams(x, p, params);

            res.succeeded = true;
        };

        // Set the external constraint function in the Optima::Problem object
        optproblem.v = [=](Optima::ConstraintResultRef res, VectorXdConstRef x, VectorXdConstRef p, VectorXdConstRef c, Optima::ConstraintOptions opts) mutable
        {
            res.val = setup.evalEquationConstraints(x, p, params);

            if(opts.eval.ddx)
                res.ddx = setup.evalEquationConstraintsGradX(x, p, params);

            if(opts.eval.ddp)
                res.ddp = setup.evalEquationConstraintsGradP(x, p, params);

            if(opts.eval.ddc)
                res.ddc = setup.evalEquationConstraintsGradParams(x, p, params);

            res.succeeded = true;
        };

        // Set the coefficient matrices Aex and Aep of the linear equality constraints
        optproblem.Aex = setup.assembleMatrixAex();
        optproblem.Aep = setup.assembleMatrixAep();

        /// Set the right-hand side vector be of the linear equality constraints.
        optproblem.be = setup.assembleVectorBe(conditions, state0);

        // Set the lower bounds of the species amounts
        optproblem.xlower = setup.assembleLowerBoundsVector(restrictions, state0);
        optproblem.xupper = setup.assembleUpperBoundsVector(restrictions, state0);
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateOptState(const ChemicalState& state0)
    {
        optstate = state0.equilibrium().optimaState();

        if(optstate.dims.x != dims.Nx)
            optstate = Optima::State(optdims);

        optstate.x.head(dims.Nn) = state0.speciesAmounts();
    }

    /// Update the chemical state object with computed optimization state.
    auto updateChemicalState(ChemicalState& state, const EquilibriumConditions& conditions)
    {
        const Params& params = conditions.params();
        const VectorXdConstRef p = optstate.p;
        const auto T = setup.extractTemperature(p, params);
        const auto P = setup.extractPressure(p, params);
        state.setTemperature(T);
        state.setPressure(P);
        state.setSpeciesAmounts(optstate.x.head(dims.Nn));
        state.equilibrium().setOptimaState(optstate);
    }

    /// Update the equilibrium sensitivity object with computed optimization sensitivity.
    auto updateEquilibriumSensitivity(EquilibriumSensitivity& sensitivity, const EquilibriumConditions& conditions)
    {
        const auto Nn = dims.Nn;
        const auto Nq = dims.Nq;
        sensitivity.dndc(optsensitivity.xc.topRows(Nn));
        sensitivity.dqdc(optsensitivity.xc.bottomRows(Nq));
        sensitivity.dpdc(optsensitivity.pc);
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

        updateOptProblem(state0, conditions, restrictions);
        updateOptState(state0);

        eqresult.optima = optsolver.solve(optproblem, optstate, optsensitivity);

        updateChemicalState(state0, conditions);
        updateEquilibriumSensitivity(sensitivity, conditions);

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
