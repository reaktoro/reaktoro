// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "Filter.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {

Filter::Filter()
{}

auto Filter::clear() -> void
{
    filter.clear();
}

auto Filter::acceptable(const Point& point) const -> bool
{
    for(const Point& x : filter)
        if(dominated(point, x))
            return false;
    return true;
}

auto Filter::extend(const Point& point) -> void
{
    // Check if the entry trying to be added to the filter is accepted to it
    if(acceptable(point)) {
        // Define the domination function to remove dominated points from the filter
        auto is_dominated = [=](const Point& x) {
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
        if(lessThan(a[i], b[i], a[i]))
            return false;
    return true;
}

} // namespace Reaktoro
