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
#include <cstddef>
#include <utility>

namespace Reaktoro {
namespace detail {

template<size_t i>
struct Counter
{
    constexpr static size_t value = i;
    constexpr operator const size_t() const { return value; }
    constexpr operator size_t() { return value; }
};

template<size_t i, size_t ibegin, size_t iend, typename Function>
constexpr auto AuxFor(Function&& f)
{
    if constexpr (i < iend) {
        f(Counter<i>{});
        AuxFor<i + 1, ibegin, iend>(std::forward<Function>(f));
    }
}

} // namespace detail

/// Generate evaluation statements `f(ibegin); f(ibegin + 1); ...; f(iend-1);` at compile time.
template<size_t ibegin, size_t iend, typename Function>
constexpr auto For(Function&& f)
{
    detail::AuxFor<ibegin, ibegin, iend>(std::forward<Function>(f));
}

/// Generate evaluation statements `f(0); f(1); ...; f(iend-1);` at compile time.
template<size_t iend, typename Function>
constexpr auto For(Function&& f)
{
    For<0, iend>(std::forward<Function>(f));
}

/// Generate evaluation statements `f(arg0); f(arg1); ...; f(argn);` at compile time.
template<typename Function>
constexpr auto ForEach(Function const& f) -> void
{
}

/// Generate evaluation statements `f(arg0); f(arg1); ...; f(argn);` at compile time.
template<typename Function, typename Arg, typename... Args>
constexpr auto ForEach(Function const& f, Arg const& arg, Args const&... args) -> void
{
    f(arg);
    ForEach(f, args...);
}

} // namespace Reaktoro
