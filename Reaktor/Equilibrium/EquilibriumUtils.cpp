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
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktor {

auto equilibrate(ChemicalState& state) -> EquilibriumResult
{
    return equilibrate(state, {});
}

auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult
{
    ChemicalSystem system = state.system();

    EquilibriumSolver solver;
    solver.setOptions(options);

    EquilibriumProblem problem(system);
    problem.setTemperature(state.temperature());
    problem.setPressure(state.pressure());
    problem.setElementAmounts(state.elementAmounts());
    problem.setCharge(0.0);

    return solver.solve(problem, state);
}

} // namespace Reaktor
