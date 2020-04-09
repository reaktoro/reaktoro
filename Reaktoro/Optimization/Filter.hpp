// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <list>
#include <vector>

namespace Reaktoro {

/// A type that describes an optimisation filter
class Filter
{
public:
    /// A type that describes an entry point of an optimisation filter
    using Point = std::vector<double>;

    /// Construc a default Filter instance
    Filter();

    /// Clear the filter by removing all its points.
    auto clear() -> void;

    /// Check if a given point is acceptable to the filter
    /// @param point The point to be checked
    auto acceptable(const Point& point) const -> bool;

    /// Extend the filter with a new point.
    /// This method not only inserts a new point to the filter,
    /// but also removes all existent dominated points with
    /// respect to the new one.
    /// @param point The point to be inserted in the filter
    auto extend(const Point& point) -> void;

    /// Check if a point is dominated by another.
    /// Check if point `a` is dominated by `b`, that is, if `a` > `b` componentwise.
    static auto dominated(const Point& a, const Point& b) -> bool;

private:
    /// The list of points in the filter
    std::list<Point> filter;
};

} // namespace Reaktoro
