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

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/// A type that describes the entry of an optimization filter
typedef std::vector<double> FilterEntry;

/// A type that describes an optimization filter
typedef std::list<FilterEntry> Filter;

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

/// Compute the largest step length @f$\alpha@f$ such that
/// @f$\mathbf{p} + \alpha\Delta\mathbf{p}@f$ is on the
/// lower bound @f$\mathbf{x}_l=\mathbf{0}@f$.
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
auto largestStep(const Vector& p, const Vector& dp) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\mathbf{p}+\alpha\Delta\mathbf{p}\geq(1-\tau)\mathbf{p}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
auto fractionToTheBoundary(const Vector& p, const Vector& dp, double tau) -> double;

/// Check if a float number is less than another by a base value.
/// The following is used for the comparison:
/// @f[a < b + 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto lessThan(double a, double b, double baseval) -> bool;

/// Check if a float number is greater than another by a base value.
/// The following is used for the comparison:
/// @f[a > b - 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto greaterThan(double a, double b, double baseval) -> bool;

} // namespace Reaktor
