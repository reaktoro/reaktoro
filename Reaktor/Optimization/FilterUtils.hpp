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

#pragma once

// C++ includes
#include <list>
#include <vector>

namespace Reaktor {

/// A type that describes the entry of an optimization filter
typedef std::vector<double> FilterEntry;

/// A type that describes an optimization filter
typedef std::list<FilterEntry> Filter;

/// Check if a filter entry is dominated by another.
/// Check if entry `a` is dominated by `b`, that is, if `a` > `b` componentwise.
/// @param a The filter entry `a`
/// @param b The filter entry `b`
auto dominated(const FilterEntry& a, const FilterEntry& b) -> bool;

/// Check if an entry is acceptable to a filter
/// @param entry The entry to be checked
/// @param filter The filter where the entry is checked
auto acceptable(const Filter& filter, const FilterEntry& entry) -> bool;

/// Extend a filter with a new entry.
/// This method not only insert a new entry to the filter,
/// but also removes all existent dominated entries with
/// the insertion of the new one.
/// @param filter The filter where the entry is inserted
/// @param entry The entry to be inserted in the filter
auto extend(Filter& filter, const FilterEntry& entry) -> void;

/// Clear a filter.
/// @param filter The filter to be cleared.
auto clear(Filter& filter) -> void;

} // namespace Reaktor
