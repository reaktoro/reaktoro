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

namespace Reaktoro {

template <typename Container, typename T>
auto index(const Container& c, const T& value) -> std::ptrdiff_t
{
    auto i = std::find(c.begin(), c.end(), value) - c.begin();
    return i < c.size() ? i : -1;
}

template <typename Container, typename Predicate>
auto indexfn(const Container& c, const Predicate& pred) -> std::ptrdiff_t
{
    auto i = std::find_if(c.begin(), c.end(), pred) - c.begin();
    return i < c.size() ? i : -1;
}

template <typename Container, typename Predicate>
auto filter(const Container& c, const Predicate& pred)
{
    Container res;
    std::copy_if(c.begin(), c.end(), std::back_inserter(res), pred);
    return res;
}

template <typename Container, typename Predicate>
auto remove(const Container& c, const Predicate& pred)
{
    return filter(c, [&](auto&& x) { return !pred(x); });
}

template <typename Container>
auto unique(const Container& c)
{
    Container res(c);
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    return res;
}

template <typename Container, typename Result, typename Function>
auto transform(const Container& c, Result& res, const Function& f)
{
    std::transform(c.begin(), c.end(), res.begin(), f);
}

template <typename Container>
auto merge(const Container& a, const Container& b)
{
    Container res(a);
    res.insert(res.end(), b.begin(), b.end());
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    return res;
}

template <typename Container, typename T>
auto contains(const Container& c, const T& value)
{
    return std::find(c.begin(), c.end(), value) != c.end();
}

template <typename Container, typename Predicate>
auto containsfn(const Container& c, const Predicate& pred)
{
    return std::find_if(c.begin(), c.end(), pred) != c.end();
}

template <typename ContainerA, typename ContainerB>
auto contained(const ContainerA& a, const ContainerB& b)
{
    for(auto const& x : a)
        if (!contains(b, x))
            return false;
    return true;
}

} // namespace Reaktoro
