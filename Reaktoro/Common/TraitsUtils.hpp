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
#include <limits>
#include <type_traits>

namespace Reaktoro {

template<bool value>
using EnableIf = std::enable_if_t<value>;

template<bool value>
using Requires = std::enable_if_t<value, bool>;

template<typename T>
using Decay = std::decay_t<T>;

template<typename T>
constexpr auto isArithmetic = std::is_arithmetic_v<T>;

template<typename T>
constexpr auto isInteger = std::numeric_limits<Decay<T>>::is_integer;

template<typename T>
constexpr auto isFloatingPoint = std::is_floating_point_v<T>;

template<typename From, typename To>
constexpr auto isConvertible = std::is_convertible_v<From, To>;

template<typename T, typename U>
constexpr auto isSame = std::is_same_v<Decay<T>, Decay<U>>;

template<typename T, typename U>
constexpr auto isBaseOf = std::is_base_of_v<Decay<T>, Decay<U>>;

namespace detail {

template<typename T>
struct isFunction { static constexpr auto value = false; };

template<typename Signature>
struct isFunction<std::function<Signature>> { static constexpr auto value = true; };

template<typename T>
struct asFunction : public asFunction<decltype(&T::operator())> {};

template<typename Ret, typename... Args>
struct asFunction<Ret(Args...)> { using type = std::function<Ret(Args...)>; };

template<typename Ret, typename... Args>
struct asFunction<Ret(*)(Args...)> { using type = std::function<Ret(Args...)>; };

template<typename Class, typename Ret, typename... Args>
struct asFunction<Ret(Class::*)(Args...) const> { using type = std::function<Ret(Args...)>; };

template<typename T, typename U, typename... Us>
constexpr auto isOneOf()
{
    if constexpr (sizeof...(Us) > 0)
        return isSame<T, U> || isOneOf<T, Us...>();
    else return isSame<T, U>;
}

} // namespace detail

template<typename T>
constexpr auto isFunction = detail::isFunction<Decay<T>>::value;

/// Convert lambda/function pointers/member functions to `std::function`.
/// This exists in Reaktoro because AppleClang 9.0/10.0/11.0 do not perform
/// template type deduction for `std::function`. Instead of just writing
/// `std::function(f)`, we need to use instead `asFunction(f)`.
template<typename Fun>
constexpr auto asFunction(const Fun& f) { return typename detail::asFunction<Fun>::type{f}; } // Reference: https://stackoverflow.com/a/39182901/418875

template<typename T, typename U, typename... Us>
constexpr auto isOneOf = detail::isOneOf<T, U, Us...>();

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
