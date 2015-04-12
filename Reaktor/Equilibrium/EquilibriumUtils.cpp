// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

// Reaktor includes
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktor {

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
    const auto& system = problem.system();
    const auto& partition = problem.partition();
    const auto be = rows(problem.elementAmounts(), partition.indicesEquilibriumElements());
    EquilibriumSolver solver(system);
    solver.setOptions(options);
    solver.setPartition(partition);
    return solver.solve(state, be);
}

} // namespace Reaktor
