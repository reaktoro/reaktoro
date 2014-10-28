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

#include "FilterUtils.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes
#include <Reaktor/Optimization/AlgorithmUtils.hpp>

namespace Reaktor {

auto dominated(const FilterEntry& a, const FilterEntry& b) -> bool
{
    for(unsigned i = 0; i < a.size(); ++i)
        if(lessThan(a[i], b[i], a[i])) return false;
    return true;
}

auto acceptable(const Filter& filter, const FilterEntry& entry) -> bool
{
    for(const FilterEntry& point : filter)
        if(dominated(entry, point)) return false;
    return true;
}

auto extend(Filter& filter, const FilterEntry& entry) -> void
{
    // Check if the entry trying to be added to the filter is accepted to it
    if(acceptable(filter, entry))
    {
        // Define the domination function to remove dominated points from the filter
        auto is_dominated = [=](const FilterEntry& point)
        {
            return dominated(point, entry);
        };

        // Remove all dominated entries in the filter
        filter.erase(std::remove_if(filter.begin(), filter.end(), is_dominated), filter.end());

        // Add the new entry to the filter
        filter.push_back(entry);
    }
}

auto clear(Filter& filter) -> void
{
    filter.clear();
}

} // namespace Reaktor
