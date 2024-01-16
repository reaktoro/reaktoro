// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <functional>
#include <utility>
#include <vector>

namespace Reaktoro {

inline void __hashCombine(std::size_t& seed) {}

template<typename T, typename... Rest>
inline void __hashCombine(std::size_t& seed, T const& v, Rest... rest)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    __hashCombine(seed, rest...);
}

/// Return the hash combine of the hash number of given values.
/// https://stackoverflow.com/a/38140932/418875
template<typename T, typename... Rest>
auto hashCombine(std::size_t seed, T const& v, Rest... rest) -> std::size_t
{
    __hashCombine(seed, v, rest...);
    return seed;
}

/// Return the hash of a vector of values.
/// https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
template<typename Vec>
auto hashVector(Vec const& vec) -> std::size_t
{
    std::size_t seed = vec.size();
    for(auto const& val : vec)
        seed = hashCombine(seed, val);
    return seed;
}

} // namespace Reaktoro

/// Specialize std::hash for `std::vector<T, A>` so that it can be used as key in `std::unordered_map`.
template<typename T, typename A>
struct std::hash<std::vector<T, A>>
{
    std::size_t operator()(std::vector<T, A> const& obj) const
    {
        return Reaktoro::hashVector(obj);
    }
};

/// Specialize std::hash for `std::pair<L, R>` so that it can be used as key in `std::unordered_map`.
template<typename L, typename R>
struct std::hash<std::pair<L, R>>
{
    std::size_t operator()(std::pair<L, R> const& obj) const
    {
        return Reaktoro::hashCombine(0, obj.first, obj.second);
    }
};
