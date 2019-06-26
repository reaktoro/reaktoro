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
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumProfiling.hpp>

namespace Reaktoro {

/// Provide profiling information of the operations in a reactive transport time step calculation.
struct ReactiveTransportProfiling
{
    /// The time spent for solving the transport calculations without chemical reactions.
    double time_transport = 0.0;

    /// The time spent for performing all chemical equilibrium calculations.
    double time_equilibrium = 0.0;

    /// The boolean flags that tell if a cell had a successful smart equilibrium calculation without learning.
    std::vector<bool> successful_smart_equilibrium_estimation_at_cell;

    /// The boolean flags that tell if a cell had a successful smart equilibrium calculation without learning.
    std::vector<EquilibriumResult> successful_smart_equilibrium_estimation_at_cell;

    /// The boolean flags that tell if a cell had a successful smart equilibrium calculation without learning.
    std::vector<SmartEquilibriumResult> successful_smart_equilibrium_estimation_at_cell;
};

} // namespace Reaktoro
