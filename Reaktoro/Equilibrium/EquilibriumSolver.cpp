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
#include <Optima/Solver.hpp>
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {
namespace detail {

/// Return the Optima::Dims object with dimensions of the optimization problem.
auto dims(const EquilibriumProblem& problem) -> Optima::Dims
{
    Optima::Dims dims;
    dims.x = problem.dims().Nx;
    dims.be = problem.dims().Nc;
    return dims;
}

/// Return a partially initialized Optima::Problem object
auto initOptProblem(const EquilibriumProblem& problem) -> Optima::Problem
{
    Optima::Problem oproblem(dims(problem));
    // TODO: This needs to be improved in Optima. Create a Constraints class
    // where constraints can be set. Implement Solver(const Constraints&).
    // Ensure the Ae matrix is not pure zeros! Use:
    //
    // solver.problem.f = new_objective_fn()
    // solver.problem.be = new_be_values()
    // solver.problem.xlower = new_lower_bounds()
    // solver.problem.Ae does not exist because cannot be changed after construction!
    // solver.solve(state)

    // Set the coefficient matrix of the linear equality constraints
    oproblem.Aex = problem.conservationMatrix();

    return oproblem;
}

} // namespace detail

struct EquilibriumSolver::Impl
{
    /// The chemical system associated with this equilibrium solver
    ChemicalSystem system;

    /// The equilibrium constraints associated with this equilibrium solver
    EquilibriumConstraints constraints;

    /// The dimensions related data in a constrained equilibrium problem
    EquilibriumDims dims;

    /// The equilibrium problem with conservation matrix and objective function
    EquilibriumProblem problem;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The optimization problem corresponding to the constrained equilibrium calculation.
    Optima::Problem optproblem;

    /// The optimization state of the calculation
    Optima::State optstate;

    /// The solver for the optimization calculations
    Optima::Solver optsolver;

    /// The auxiliary vector to store the amounts of the species
    ArrayXd n0;

    /// Construct a Impl instance with given EquilibriumConstraints object
    Impl(const EquilibriumConstraints& ecs)
    : system(ecs.system()),
      constraints(ecs),
      dims(ecs),
      problem(ecs),
      optproblem(detail::initOptProblem(problem)),
      optstate(detail::dims(problem)),
      optsolver(optproblem)
    {
        // Prevent new control variables, functional constraints, and inert reactions
        constraints.lock();

        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    /// Set the options of the equilibrium solver.
    auto setOptions(const EquilibriumOptions& options_) -> void
    {
        // Update the options of the equilibrium calculation
        options = options_;

        // Ensure some options have proper values
        error(options.epsilon <= 0, "EquilibriumOptions::epsilon cannot be zero or negative.");

        // Set the lower bounds of the species amounts
        optproblem.xlower.head(dims.Nn).fill(options.epsilon);

        // Initialize the names of the primal and dual variables
        if(options.optima.output.active)
        {
            // Use `n` instead of `x` to name the variables
            options.optima.output.xprefix = "n";

            // Define some auxiliary references to the variables names
            auto& xnames = options.optima.output.xnames;

            // Initialize the names of the primal variables `n`
            for(auto species : system.species())
                xnames.push_back(species.name());
        }

        // Set the options of the optimization solver
        optsolver.setOptions(options.optima);
    }

    /// Update the optimization problem before a new equilibrium calculation.
    auto updateOptProblem(ChemicalState& state0, ArrayXdConstRef b)
    {
        // Update the equilibrium problem with updated equilibrium constraints
        problem.update(constraints);

        // Correct zero amounts of species before calculating properties
        n0 = state0.speciesAmounts();
        n0 = (n0 <= 0.0).select(options.epsilon, n0);
        state0.setSpeciesAmounts(n0);

        // Ensure the chemical properties in state0 are in sync with its (T,P,n)
        state0.props().update(state0); // TODO: This needs to be improved to avoid accidental unsync issues!

        // Create the objective function corresponding to the constrained equilibrium calculation
        EquilibriumObjective obj = problem.objective(state0);

        // Update the objective function in the Optima::Problem object
        optproblem.f = [=](VectorXdConstRef x, VectorXdConstRef p, Optima::ObjectiveResult& res)
        {
            if(res.requires.f) res.f = obj.f(x);
            if(res.requires.fx) obj.g(x, res.fx);
            if(res.requires.fxx) obj.H(x, res.fxx);
        };

        /// Update the right-hand side vector of the linear equality constraints
        optproblem.be = b;

        // Update the lower and upper bounds of the variables
        problem.xlower(state0, optproblem.xlower);
        problem.xupper(state0, optproblem.xupper);
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateOptState(const ChemicalState& state0)
    {
        const auto& n = state0.speciesAmounts();
        const auto& y = state0.equilibrium().lagrangeMultipliers();
        const auto& s = state0.equilibrium().complementarityVariables();
        const auto& v = state0.equilibrium().controlVariables();

        optstate.x.head(dims.Nn) = n;

        if(v.size() == dims.Ncv) optstate.x.tail(dims.Ncv) = v;
        if(y.size() == dims.Nc)  optstate.y = y;
        if(s.size() == dims.Ne)  optstate.s.head(dims.Nn) = s;

        optstate.s.tail(dims.Ncv).fill(0.0); // because the control variables don't have lower/upper bounds
    }

    /// Update the initial state variables before the new equilibrium calculation.
    auto updateChemicalState(ChemicalState& state)
    {
        state.speciesAmounts() = optstate.x.head(dims.Nn);
        state.equilibrium().controlVariables() = optstate.x.tail(dims.Ncv);
        state.equilibrium().lagrangeMultipliers() = optstate.y;
        state.equilibrium().complementarityVariables() = optstate.s.head(dims.Nn);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0) -> EquilibriumResult
    {
        // The conservation matrix Cn in C = [Cn Cp Cq] corresponding to the
        // chemical species (not the control variables!)
        const auto Cn = optproblem.Aex.leftCols(dims.Nn);

        // The initial amounts of the species
        n0 = state0.speciesAmounts();

        // Compute the amounts of the components that need to be conserved
        optproblem.be = Cn * n0.matrix();

        // Solve the equilibrium problem
        return solve(state0, optproblem.be);
    }

    /// Solve an equilibrium problem with given chemical state in disequilibrium.
    auto solve(ChemicalState& state0, ArrayXdConstRef b) -> EquilibriumResult
    {
        EquilibriumResult eqresult;

        updateOptProblem(state0, b);
        updateOptState(state0);

        eqresult.optima = optsolver.solve(optstate, optproblem);

        updateChemicalState(state0);

        return eqresult;
    }
};

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(EquilibriumConstraints(system)))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumConstraints& constraints)
: pimpl(new Impl(constraints))
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

auto EquilibriumSolver::constraints() const -> const EquilibriumConstraints&
{
    return pimpl->constraints;
}

auto EquilibriumSolver::constraints() -> EquilibriumConstraints&
{
    return pimpl->constraints;
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return pimpl->solve(state);
}

} // namespace Reaktoro
