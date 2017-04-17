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

namespace Reaktoro {

template<typename T>
auto index(const T& value, const std::vector<T>& values) -> Index
{
    return std::find(values.begin(), values.end(), value) - values.begin();
}

inline auto index(const std::string& word, const std::vector<std::string>& strings) -> Index
{
    return index<std::string>(word, strings);
}

template<typename NamedValues>
auto index(const std::string& name, const NamedValues& values) -> Index
{
    Index idx = 0;
    for(const auto& value : values)
        if(value.name() == name) return idx; else ++idx;
    return idx;
}

template<typename NamedValue, typename NamedValues>
auto index(const NamedValue& value, const NamedValues& values) -> Index
{
    return index(value.name(), values);
}

template<typename Names, typename NamedValues>
auto indexAny(const Names& names, const NamedValues& values) -> Index
{
    for(auto& name : names)
    {
        const Index i = index(name, values);
        if(i < values.size())
            return i;
    }
    return values.size();
}

template<typename NamedValues>
auto indices(const std::vector<std::string>& names, const NamedValues& values) -> Indices
{
    Indices idxs;
    idxs.reserve(names.size());
    for(const auto& name : names)
        idxs.push_back(index(name, values));
    return idxs;
}

template<typename NamedValues>
auto indices(const NamedValues& subvalues, const NamedValues& values) -> Indices
{
    Indices idxs;
    idxs.reserve(subvalues.size());
    for(const auto& value : subvalues)
        idxs.push_back(index(value, values));
    return idxs;
}

inline auto indices(const std::vector<std::string>& words, const std::vector<std::string>& strings) -> Indices
{
    Indices indices;
    indices.reserve(words.size());
    for(const std::string iter : words)
        indices.push_back(index(iter, strings));

    return indices;
}

template<typename NamedValues>
auto contained(const std::string& name, const NamedValues& values) -> bool
{
    return index(name, values) < values.size();
}

template<typename Container, typename = typename std::enable_if<!std::is_same<typename Container::value_type, std::string>::value>::type>
auto contained(const typename Container::value_type& value, const Container& values) -> bool
{
    return std::count(values.begin(), values.end(), value);
}

template<typename Container>
auto contained(const Container& values1, const Container& values2) -> bool
{
    for(const auto& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

template<typename T>
auto unify(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::set<T> set;

    set.insert(values1.begin(), values1.end());
    set.insert(values2.begin(), values2.end());

    return std::vector<T>(set.begin(), set.end());
}

template<typename T>
auto intersect(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::vector<T> intersection;

    for(const T& value : values1)
        if(contained(value, values2))
            intersection.push_back(value);

    return intersection;
}

template<typename T>
auto difference(const std::vector<T>& values1, const std::vector<T>& values2) -> std::vector<T>
{
    std::vector<T> diff;

    for(const T& value : values1)
        if(!contained(value, values2))
            diff.push_back(value);

    return diff;
}

template<typename T>
auto emptyIntersection(const std::vector<T>& values1, const std::vector<T>& values2) -> bool
{
    for(const T& value : values1)
        if(contained(value, values2))
            return false;
    return true;
}

template<typename T>
auto emptyDifference(const std::vector<T>& values1, const std::vector<T>& values2) -> bool
{
    for(const T& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

template<typename Container>
auto equal(const Container& values1, const Container& values2) -> bool
{
    if(values1.size() != values2.size())
        return false;
    for(const auto& value : values1)
        if(!contained(value, values2))
            return false;
    return true;
}

template<typename Container>
auto isunique(Container values) -> bool
{
    std::set<Index> tmp(values.begin(), values.end());
    return tmp.size() == values.size();
}

template<typename T>
auto unique(std::vector<T> values) -> std::vector<T>
{
    auto it = std::unique(values.begin(), values.end());

    values.resize(std::distance(values.begin(), it));

    return values;
}

template<typename T>
auto range(T first, T last, T step) -> std::vector<T>
{
    unsigned size = unsigned((last - first)/step);
    std::vector<T> range(size);
    for(unsigned i = 0; i < size; ++i)
        range[i] = first + i*step;
    return range;
}

template<typename T>
auto range(T first, T last) -> std::vector<T>
{
    return range(first, last, static_cast<T>(1));
}

template<typename T>
auto range(T last) -> std::vector<T>
{
    return range(static_cast<T>(0), last, static_cast<T>(1));
}

template<typename T, typename Predicate>
auto filter(const std::vector<T>& values, Predicate predicate) -> std::vector<T>
{
    std::vector<T> list;
    for(const T& value : values)
        if(predicate(value))
            list.push_back(value);
    return list;
}

template<typename T, typename Predicate>
auto remove(const std::vector<T>& values, Predicate predicate) -> std::vector<T>
{
    std::vector<T> list;
    for(const T& value : values)
        if(!predicate(value))
            list.push_back(value);
    return list;
}

template<typename T>
auto extract(const std::vector<T>& values, const Indices& indices) -> std::vector<T>
{
    std::vector<T> extracted_values;

    for(const auto& idx : indices)
        extracted_values.push_back(values[idx]);

    return extracted_values;
}

template<typename Container>
auto min(const Container& values) -> typename Container::value_type
{
    return *std::min_element(values.begin(), values.end());
}

template<typename Container>
auto max(const Container& values) -> typename Container::value_type
{
    return *std::max_element(values.begin(), values.end());
}

} // namespace Reaktoro


