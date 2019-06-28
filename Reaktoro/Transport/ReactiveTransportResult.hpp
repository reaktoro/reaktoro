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

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>

namespace Reaktoro {

/// Provide timing information of the operations in a reactive transport time step calculation.
struct ReactiveTransportTiming
{
    /// The time spent during the last time step of a reactive transport calculation.
    double step = 0.0;

    /// The time spent during the last time step for solving the transport equations, without chemical reactions.
    double transport = 0.0;

    /// The time spent during the last time step for performing chemical equilibrium calculations.
    double equilibrium = 0.0;
};

/// Provide result information of a reactive transport time step calculation.
struct ReactiveTransportResult
{
    /// The result of each fluid element transport calculation during the last time step.
    std::vector<TransportResult> transport_of_element;

    /// The result of the equilibrium calculation in each cell during the last time step.
    std::vector<EquilibriumResult> equilibrium_at_cell;

    /// The result of the smart equilibrium calculation in each cell during the last time step.
    std::vector<SmartEquilibriumResult> smart_equilibrium_at_cell;

    /// The timing information of the operations in a reactive transport time step calculation.
    ReactiveTransportTiming timing;
};

} // namespace Reaktoro
