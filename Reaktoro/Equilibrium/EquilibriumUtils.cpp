// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
#include "EquilibriumUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {
namespace {

auto equilibrateAux(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumProblem& problem, EquilibriumOptions options) -> EquilibriumResult
{
    // Define auxiliary references to problem data
    const auto& system = problem.system();
    const auto& partition = problem.partition();
    const auto& iee = partition.indicesEquilibriumElements();
    const auto& T = problem.temperature();
    const auto& P = problem.pressure();
    const auto& b = problem.elementAmounts();
    const auto& be = rows(b, iee);

    // Set the temperature and pressure of the chemical state
    state.setTemperature(T);
    state.setPressure(P);

    // The result of the equilibrium calculation
    EquilibriumResult res;

    // Check if the equilibrium problem is in fact an inverse equilibrium problem
    if(problem.isInverseProblem())
    {
        // Perform an inverse equilibrium calculation
        EquilibriumInverseSolver solver(system);
        solver.setPartition(partition);
        solver.setOptions(options);
        res = solver.solve(state, problem);
        sensitivity = solver.sensitivity();
    }
    else
    {
        // Perform a direct equilibrium calculation
        EquilibriumSolver solver(system);
        solver.setPartition(partition);
        solver.setOptions(options);
        res = solver.solve(state, T, P, be);
        sensitivity = solver.sensitivity();
    }

    Assert(res.optimum.succeeded, "Could not calculate the equilibrium state of the system.",
        "Convergence could not be established with given equilibrium conditions, initial guess, and/or numerical parameters.");

    return res;
}

} // namespace

auto equilibrate(ChemicalState& state) -> EquilibriumResult
{
    EquilibriumOptions options;
    return equilibrate(state, options);
}

auto equilibrate(ChemicalState& state, const Partition& partition) -> EquilibriumResult
{
    return equilibrate(state, partition, {});
}

auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult
{
    ChemicalSystem system = state.system();

    return equilibrate(state, Partition(system), options);
}

auto equilibrate(ChemicalState& state, const Partition& partition, const EquilibriumOptions& options) -> EquilibriumResult
{
    ChemicalSystem system = state.system();

    EquilibriumProblem problem(system);
    problem.setPartition(partition);
    problem.setTemperature(state.temperature());
    problem.setPressure(state.pressure());
    problem.setElementAmounts(state.elementAmounts());

    return equilibrate(state, problem, options);
}

auto equilibrate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return equilibrate(state, problem, {});
}

auto equilibrate(ChemicalState& state, const EquilibriumProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult
{
    EquilibriumSensitivity sensitivity;
    return equilibrateAux(state, sensitivity, problem, options);
}

auto equilibrate(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return equilibrate(state, sensitivity, problem, {});
}

auto equilibrate(ChemicalState& state, EquilibriumSensitivity& sensitivity, const EquilibriumProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult
{
    return equilibrateAux(state, sensitivity, problem, options);
}

auto equilibrate(const EquilibriumProblem& problem) -> ChemicalState
{
    return equilibrate(problem, {});
}

auto equilibrate(const EquilibriumProblem& problem, const EquilibriumOptions& options) -> ChemicalState
{
    ChemicalState state(problem.system());
    equilibrate(state, problem, options);
    return state;
}

} // namespace Reaktoro
