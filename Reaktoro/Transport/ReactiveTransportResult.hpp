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

namespace Reaktoro {

struct Results
{
    /// Total CPU time required by smart equilibrium scheme
    double smart_total;

    /// Total CPU time required by smart equilibrium scheme
    /// excluding the costs for the search of the closest reference states.
    double smart_total_ideal_search;

    /// Total CPU time required by smart equilibrium scheme
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_total_ideal_search_store;

    /// Total CPU time required by conventional equilibrium scheme
    double conv_total;
};

/// Use this class to collect modeling results per one step of reactive transport.
struct ReactiveTransportResult
{
    /// The boolean flag that indicates if the reactive transport time step was successful.
    bool successful = false;

    /// Flag for the smart equilibrium
    bool smart;

    /// Contains reactive transport tracking results
    EquilibriumResult equilibrium;

    /// The indices of the cells where smart equilibrium calculation was successful
    std::vector<Index> equilibrium_smart_successfull_cells;
};

} // namespace Reaktoro
