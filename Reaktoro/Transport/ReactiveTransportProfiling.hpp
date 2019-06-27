// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#pragma once

// C++ includes
#include <vector>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumProfiling.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumProfiling.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {

/// Provide profiling information of the operations in a reactive transport time step calculation.
struct ReactiveTransportProfiling
{
    /// The time spent for solving the transport calculations without chemical reactions.
    double time_transport = 0.0;

    /// The time spent for performing all chemical equilibrium calculations.
    double time_equilibrium = 0.0;

    /// The result of each cell equilibrium calculation using EquilibriumSolver in a time step.
    std::vector<EquilibriumResult> equilibrium_result_at_cell;

    /// The result of each cell equilibrium calculation using SmartEquilibriumSolver in a time step.
    std::vector<SmartEquilibriumResult> smart_equilibrium_result_at_cell;

    /// The result of each cell equilibrium calculation using EquilibriumSolver in a time step.
    std::vector<EquilibriumProfiling> equilibrium_profiling_at_cell;

    /// The result of each cell equilibrium calculation using SmartEquilibriumSolver in a time step.
    std::vector<SmartEquilibriumProfiling> smart_equilibrium_profiling_at_cell;
};

} // namespace Reaktoro
