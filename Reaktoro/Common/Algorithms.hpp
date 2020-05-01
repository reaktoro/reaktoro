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
#include <algorithm>
#include <cassert>
#include <vector>

namespace Reaktoro {

/// @def RKT_LAMBDA(x, expr)
/// A convenient macro that expands to `[&](auto&& x) { return expr; }`.
#define RKT_LAMBDA(x, expr) [&](const auto& x) { return expr; }

/// Return the index of item `x` in container `c` or the number of items if not found.
template<typename Container, typename T>
auto index(const Container& c, const T& x) -> std::size_t
{
    return std::find(c.begin(), c.end(), x) - c.begin();
}

/// Return the index of item `x` in container `c` for which `pred(x)` evaluates to true or the number of items if not found.
template<typename Container, typename Predicate>
auto indexfn(const Container& c, const Predicate& pred) -> std::size_t
{
    return std::find_if(c.begin(), c.end(), pred) - c.begin();
}

/// Return a container with items `x` for which `pred(x)` evaluates to true.
template<typename Container, typename Predicate>
auto filter(const Container& c, const Predicate& pred)
{
    Container res;
    std::copy_if(c.begin(), c.end(), std::back_inserter(res), pred);
    return res;
}

/// Return a container without items `x` for which `pred(x)` evaluates to true.
template<typename Container, typename Predicate>
auto remove(const Container& c, const Predicate& pred)
{
    return filter(c, [&](auto&& x) { return !pred(x); });
}

/// Return a container without duplicates.
template<typename Container>
auto unique(const Container& c)
{
    Container res(c);
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    return res;
}

/// Apply a function `f` on every item in container `c` and store in `res`.
template<typename Container, typename Result, typename Function>
auto transform(const Container& c, Result& res, const Function& f)
{
    assert(c.size() == res.size());
    std::transform(c.begin(), c.end(), res.begin(), f);
}

/// Return the items in a container at given indices.
template<typename Container, typename Indices>
auto extract(const Container& a, const Indices& indices)
{
    Container res;
    res.reserve(indices.size());
    for(auto i : indices)
        res.push_back(a[i]);
    return res;
}

/// Return a vector by applying function `f` on every item in container `c`.
template<typename Container, typename Function>
auto vectorize(const Container& c, const Function& f)
{
    using X = typename Container::value_type;
    using T = std::invoke_result_t<Function, X>;
    std::vector<T> res;
    res.resize(c.size());
    transform(c, res, f);
    return res;
}

/// Return true if container `a` contains item `x`.
template<typename Container, typename T>
auto contains(const Container& c, const T& x)
{
    return std::find(c.begin(), c.end(), x) != c.end();
}

/// Return true if container `a` contains item `x` for which `pred(x)` evaluates to true.
template<typename Container, typename Predicate>
auto containsfn(const Container& c, const Predicate& pred)
{
    return std::find_if(c.begin(), c.end(), pred) != c.end();
}

/// Return true if items in container `a` are also in container `b`.
template<typename ContainerA, typename ContainerB>
auto contained(const ContainerA& a, const ContainerB& b)
{
    for(auto const& x : a)
        if(!contains(b, x))
            return false;
    return true;
}

/// Return a container with items from both `a` and `b`.
template<typename Container>
auto concatenate(const Container& a, const Container& b)
{
    Container res(a);
    res.insert(res.end(), b.begin(), b.end());
    return res;
}

/// Return a container with items from both `a` and `b` without duplicates.
template<typename Container>
auto merge(const Container& a, const Container& b)
{
    Container res(a);
    res.insert(res.end(), b.begin(), b.end());
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    return res;
}

/// Return the intersection of two containers.
template<typename Container>
auto intersect(const Container& a, const Container& b)
{
    return filter(a, RKT_LAMBDA(x, contains(b, x)));
}

/// Return the difference of two containers.
template<typename Container>
auto difference(const Container& a, const Container& b)
{
    return filter(a, RKT_LAMBDA(x, !contains(b, x)));
}

/// Return true if containers `a` and `b` have distinct items.
template<typename ContainerA, typename ContainerB>
auto disjoint(const ContainerA& a, const ContainerB& b)
{
    for(auto&& x : a)
        if(contains(b, x))
            return false;
    return true;
}

/// Return true if containers `a` and `b` have identical items.
template<typename ContainerA, typename ContainerB>
auto identical(const ContainerA& a, const ContainerB& b)
{
    return contained(a, b) && contained(b, a);
}

/// Return a vector with given range of values and step between them.
template<typename T>
auto range(T first, T last, T step)
{
    if(last <= first) return std::vector<T>{};
    auto size = std::size_t((last - first - 1)/step) + 1;
    std::vector<T> res(size);
    for(auto i = 0; i < size; ++i)
        res[i] = first + i*step;
    return res;
}

/// Return a vector with given range of values with unit step.
template<typename T>
auto range(T first, T last)
{
    return range(first, last, static_cast<T>(1));
}

/// Return a vector with given range of values with unit step and starting from 0.
template<typename T>
auto range(T last)
{
    return range(static_cast<T>(0), last, static_cast<T>(1));
}

/// Return true if `x` is equal to at least one of the other arguments.
template<typename X, typename X0, typename... XS>
auto oneof(const X x, const X0& x0, const XS&... xs)
{
    if constexpr(sizeof...(XS) > 0)
        return x == x0 || oneof(x, xs...);
    return x == x0;
}

} // namespace Reaktoro
