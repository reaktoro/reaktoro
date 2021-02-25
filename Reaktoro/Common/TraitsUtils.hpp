// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <type_traits>

namespace Reaktoro {

namespace detail {

template<typename T>
struct isFunction { static constexpr auto value = false; };

template<typename Signature>
struct isFunction<std::function<Signature>> { static constexpr auto value = true; };

} // namespace detail

template<bool value>
using EnableIf = std::enable_if_t<value>;

template<typename T>
using Decay = std::decay_t<T>;

template<typename T>
constexpr auto isArithmetic = std::is_arithmetic_v<Decay<T>>;

template<typename T, typename U>
constexpr auto isSame = std::is_same_v<Decay<T>, Decay<U>>;

namespace detail {

template<typename T, typename U, typename... Us>
constexpr auto isOneOf()
{
    if constexpr (sizeof...(Us))
        return isSame<T, U> || isOneOf<T, Us...>();
    else return isSame<T, U>;
}

} // namespace detail

template<typename T, typename U, typename... Us>
constexpr auto isOneOf = detail::isOneOf<T, U, Us...>();

template<typename T>
constexpr auto isFunction = detail::isFunction<Decay<T>>::value;

//======================================================================
// Reference type traits
//======================================================================
namespace detail { template<typename T> struct Ref { using type = T&; }; }

template<typename T>
using Ref = typename detail::Ref<T>::type;

#define REAKTORO_DEFINE_REFERENCE_TYPE_OF(basetype, reftype) \
    namespace detail { template<> struct Ref<basetype> { using type = reftype; }; }

//======================================================================

} // namespace Reaktoro
