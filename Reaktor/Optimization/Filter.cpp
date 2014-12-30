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

#include <Reaktor/Optimization/Utils.hpp>
#include "Filter.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes

namespace Reaktor {

Filter::Filter()
{}

auto Filter::clear() -> void
{
    filter.clear();
}

auto Filter::acceptable(const Point& point) const -> bool
{
    for(const Point& x : filter)
        if(dominated(point, x)) return false;
    return true;
}

auto Filter::extend(const Point& point) -> void
{
    // Check if the entry trying to be added to the filter is accepted to it
    if(acceptable(point))
    {
        // Define the domination function to remove dominated points from the filter
        auto is_dominated = [=](const Point& x)
        {
            return dominated(x, point);
        };

        // Remove all dominated entries in the filter
        filter.erase(std::remove_if(filter.begin(), filter.end(), is_dominated), filter.end());

        // Add the new entry to the filter
        filter.push_back(point);
    }
}

auto Filter::dominated(const Point& a, const Point& b) -> bool
{
    for(unsigned i = 0; i < a.size(); ++i)
        if(lessThan(a[i], b[i], a[i])) return false;
    return true;
}

} // namespace Reaktor
