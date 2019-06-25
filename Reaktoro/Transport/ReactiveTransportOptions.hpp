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
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

/// The options for the reactive transport calculations.
struct ReactiveTransportOptions
{
    /// The options for the equilibrium calculations.
    EquilibriumOptions equilibrium;

    /// The boolean flag that indicates whether SmartEquilibriumSolver should be used.
    bool use_smart_equilibrium_solver = false;
};

} // namespace Reaktoro
