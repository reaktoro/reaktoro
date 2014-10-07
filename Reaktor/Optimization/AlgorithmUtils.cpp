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

#include "AlgorithmUtils.hpp"

// C++ includes
#include <algorithm>

namespace Reaktor {
namespace {

bool operator<(const FilterEntry& a, const FilterEntry& b)
{
    for(unsigned i = 0; i < a.size(); ++i)
        if(b[i] < a[i]) return false;
    return true;
}

bool operator>(const FilterEntry& a, const FilterEntry& b)
{
    return b < a;
}

} // namespace

auto acceptable(const FilterEntry& entry, const Filter& filter) -> bool
{
    for(const FilterEntry& point : filter)
        if(entry > point)
            return false;
    return true;
}

auto extend(const FilterEntry& entry, Filter& filter) -> void
{
    // Define the domination function to remove dominated points from the filter
    auto dominated = [=](const FilterEntry& point)
    {
        return entry < point;
    };

    // Remove all dominated entries in the filter
    filter.erase(std::remove_if(filter.begin(), filter.end(), dominated), filter.end());

    // Add the new entry to the filter
    filter.push_back(entry);
}

auto largestStep(const Vector& p, const Vector& dp) -> double
{
    Vector res = -p/dp;
    double alpha = INFINITY;
    for(unsigned i = 0; i < res.size(); ++i)
        if(res[i] > 0.0 and res[i] < alpha)
            alpha = res[i];
    return alpha;
}

auto lessThan(double lhs, double rhs, double baseval) -> bool
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    return lhs < rhs + 10.0 * epsilon * std::abs(baseval);
}

auto greaterThan(double lhs, double rhs, double baseval) -> bool
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    return lhs > rhs - 10.0 * epsilon * std::abs(baseval);
}

} // namespace Reaktor
