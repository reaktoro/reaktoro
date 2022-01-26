// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "EquilibriumUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {

auto equilibrate(ChemicalState& state) -> EquilibriumResult
{
    EquilibriumOptions options;
    return equilibrate(state, options);
}

auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult
{
    EquilibriumSolver solver(state.system());
    solver.setOptions(options);
    return solver.solve(state);
}

auto equilibrate(Material& material) -> Tuple<ChemicalState, EquilibriumResult>
{
    EquilibriumOptions options;
    return equilibrate(material, options);
}

auto equilibrate(Material& material, const EquilibriumOptions& options) -> Tuple<ChemicalState, EquilibriumResult>
{
    ChemicalState state(material.system());

    EquilibriumOptions opts(options);

    EquilibriumSolver solver(state.system());

    opts.use_ideal_activity_models = true; // force ideal activity models for the first computation
    solver.setOptions(opts);

    auto result = solver.solve(state);

    opts.use_ideal_activity_models = options.use_ideal_activity_models; // for the second computation, use what user wants (maybe ideal model again, in which case the calculation below will converge immediately).
    solver.setOptions(opts);

    result += solver.solve(state);

    return {state, result};
}

} // namespace Reaktoro
