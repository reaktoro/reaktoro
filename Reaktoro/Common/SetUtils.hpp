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

#pragma once

// C++ includes
#include <algorithm>
#include <set>
#include <vector>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>

namespace Reaktoro {

/// Find the index of a value in a container of values
template<typename T>
auto index(const T& value, const std::vector<T>& values) -> Index;

/// Find the index of the @c word in the container of @c strings
auto index(const std::string& word, const std::vector<std::string>& strings) -> Index;

/// Return the index of an entry in a container.
template<typename NamedValues>
auto index(const std::string& name, const NamedValues& values) -> Index;

/// Return the index of an entry in a container.
template<typename NamedValue, typename NamedValues>
auto index(const NamedValue& value, const NamedValues& values) -> Index;

/// Return the index of the first entry in a container of named values with any of the given names.
template<typename Names, typename NamedValues>
auto indexAny(const Names& names, const NamedValues& values) -> Index;

/// Return the indices of some entries in a container.
template<typename NamedValues>
auto indices(const std::vector<std::string>& names, const NamedValues& values) -> Indices;

/// Return the indices of some entries in a container.
template<typename NamedValues>
auto indices(const NamedValues& subvalues, const NamedValues& values) -> Indices;

/// Find the indices of the `words` in the container of `strings`
auto indices(const std::vector<std::string>& words, const std::vector<std::string>& strings) -> Indices;

/// Check if a value is contained in a container of values
template<typename Container>
auto contained(const typename Container::value_type& value, const Container& values) -> bool;

/// Check if a container of values is contained in another
template<typename Container>
auto contained(const Container& values1, const Container& values2) -> bool;

/// Return true if a named value is in a set of values.
template<typename NamedValues>
auto contained(const std::string& name, const NamedValues& values) -> bool;

/// Determine the union of two containers
template<typename T>
auto unify(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>;

/// Determine the intersection of two containers
template<typename T>
auto intersect(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>;

/// Determine the difference of two containers
template<typename T>
auto difference(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>;

/// Check if two containers have empty intersection
template<typename T>
auto emptyIntersection(const std::vector<T>& values1, const std::vector<T>& values2) -> bool;

/// Check if two containers have empty difference
template<typename T>
auto emptyDifference(const std::vector<T>& values1, const std::vector<T>& values2) -> bool;

/// Check if two containers are equal
template<typename Container>
auto equal(const Container& values1, const Container& values2) -> bool;

/// Check if a container has unique values
template<typename Container>
auto isunique(Container values) -> bool;

/// Create a container with unique values from another
template<typename T>
auto unique(std::vector<T> values) -> std::vector<T>;

/// Return a range of values
/// @param begin The begin of the sequence
/// @param end The past-the-end entry of the sequence
/// @param step The step of the sequence
template<typename T>
auto range(T first, T last, T step) -> std::vector<T>;

/// Return a range of values with unit step
/// @param begin The begin of the sequence
/// @param end The past-the-end entry of the sequence
template<typename T>
auto range(T first, T last) -> std::vector<T>;

/// Return a range of values starting from zero and increasing by one
/// @param The size of the sequence
template<typename T>
auto range(T last) -> std::vector<T>;

/// Filter the values that pass on the predicate
template<typename T, typename Predicate>
auto filter(const std::vector<T>& values, Predicate predicate) -> std::vector<T>;

/// Remove the values that pass on the predicate
template<typename T, typename Predicate>
auto remove(const std::vector<T>& values, Predicate predicate) -> std::vector<T>;

/// Extract values from a vector given a list of indices
/// @param values The values from which a extraction will be performed
/// @param indices The indices of the values to be extracted
/// @return The extracted values
template<typename T>
auto extract(const std::vector<T>& values, const Indices& indices) -> std::vector<T>;

/// Return the minimum value in a STL compatible container.
template<typename Container>
auto minValue(const Container& values) -> typename Container::value_type;

/// Return the maximum value in a STL compatible container.
template<typename Container>
auto maxValue(const Container& values) -> typename Container::value_type;

} // namespace Reaktoro

#include <Reaktoro/Common/SetUtils.hxx>
