// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See The
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <algorithm>
#include <set>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

/// Find the index of a value in a container of values
template<typename T>
inline auto index(const T& value, const std::vector<T>& values) -> Index
{
    return std::find(values.begin(), values.end(), value) - values.begin();
}

/// Find the index of the @c word in the container of @c strings
inline auto index(const std::string& word, const std::vector<std::string>& strings) -> Index
{
    return index<std::string>(word, strings);
}

/// Find the indices of the @c words in the container of @c strings
inline auto indices(const std::vector<std::string>& words, const std::vector<std::string>& strings) -> Indices
{
    Indices indices;
    indices.reserve(words.size());
    for(const std::string iter : words)
        indices.push_back(index(iter, strings));

    return indices;
}

/// Check if a value is contained in a container of values
template<typename Container>
inline auto contained(const typename Container::value_type& value, const Container& values) -> bool
{
    return std::count(values.begin(), values.end(), value);
}

/// Check if a container of values is contained in another
template<typename Container>
inline auto contained(const Container& values1, const Container& values2) -> bool
{
    for(const auto& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

/// Determine the union of two containers
template<typename T>
inline auto unify(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::set<T> set;

    set.insert(values1.begin(), values1.end());
    set.insert(values2.begin(), values2.end());

    return std::vector<T>(set.begin(), set.end());
}

/// Determine the intersection of two containers
template<typename T>
inline auto intersect(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::vector<T> intersection;

    for(const T& value : values1)
        if(contained(value, values2))
            intersection.push_back(value);

    return intersection;
}

/// Determine the difference of two containers
template<typename T>
inline auto difference(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::vector<T> diff;

    for(const T& value : values1)
        if(!contained(value, values2))
            diff.push_back(value);

    return diff;
}

/// Check if two containers have empty intersection
template<typename T>
inline auto emptyIntersection(const std::vector<T>& values1, const std::vector<T>& values2) -> bool
{
    for(const T& value : values1)
        if(contained(value, values2))
            return false;
    return true;
}

/// Check if two containers have empty difference
template<typename T>
inline auto emptyDifference(const std::vector<T>& values1, const std::vector<T>& values2) -> bool
{
    for(const T& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

/// Check if two containers are equal
template<typename Container>
inline auto equal(const Container& values1, const Container& values2) -> bool
{
    if(values1.size() != values2.size())
        return false;
    for(const auto& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

/// Check if a container has unique values
template<typename Container>
inline auto isUnique(Container values) -> bool
{
    std::set<Index> tmp(values.begin(), values.end());
    return tmp.size() == values.size();
}

/// Create a container with unique values from another
template<typename T>
inline auto unique(std::vector<T> values) -> std::vector<T>
{
    auto it = std::unique(values.begin(), values.end());

    values.resize(std::distance(values.begin(), it));

    return values;
}

/// Return a range of values
/// @param begin The begin of the sequence
/// @param end The past-the-end entry of the sequence
/// @param step The step of the sequence
template<typename T>
inline auto range(T first, T last, T step) -> std::vector<T>
{
    unsigned size = unsigned((last - first)/step);
    std::vector<T> range(size);
    for(unsigned i = 0; i < size; ++i)
        range[i] = first + i*step;
    return range;
}

/// Return a range of values with unit step
/// @param begin The begin of the sequence
/// @param end The past-the-end entry of the sequence
template<typename T>
inline auto range(T first, T last) -> std::vector<T>
{
    return range(first, last, static_cast<T>(1));
}

/// Return a range of values starting from zero and increasing by one
/// @param The size of the sequence
template<typename T>
inline auto range(T last) -> std::vector<T>
{
    return range(static_cast<T>(0), last, static_cast<T>(1));
}

/// Filter the values that pass on the predicate
template<typename T, typename Predicate>
inline auto filter(const std::vector<T>& values, Predicate predicate) -> std::vector<T>
{
    std::vector<T> filtered_values;

    for(const T& value : values)
        if(predicate(value))
            filtered_values.push_back(value);

    return filtered_values;
}

/// Remove the values that pass on the predicate
template<typename T, typename Predicate>
inline auto remove(const std::vector<T>& values, Predicate predicate) -> std::vector<T>
{
    std::vector<T> filtered_values;

    for(const T& value : values)
        if(!predicate(value))
            filtered_values.push_back(value);

    return filtered_values;
}

/// Extract values from a vector given a list of indices
/// @param values The values from which a extraction will be performed
/// @param indices The indices of the values to be extracted
/// @return The extracted values
template<typename T>
inline auto extract(const std::vector<T>& values, const Indices& indices) -> std::vector<T>
{
    std::vector<T> extracted_values;

    for(const auto& idx : indices)
        extracted_values.push_back(values[idx]);

    return extracted_values;
}

/// Return the minimum value in a STL compatible container.
template<typename Container>
auto min(const Container& values) -> typename Container::value_type
{
    return *std::min_element(values.begin(), values.end());
}

/// Return the maximum value in a STL compatible container.
template<typename Container>
auto max(const Container& values) -> typename Container::value_type
{
    return *std::max_element(values.begin(), values.end());
}

} // namespace Reaktoro


