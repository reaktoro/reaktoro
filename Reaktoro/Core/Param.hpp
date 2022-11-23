// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <limits>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NumberTraits.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// A type used to represent the value of a parameter and its lower and upper bounds.
class Param
{
public:
    /// Construct a default Param object.
    Param();

    /// Construct a Param object with given value.
    Param(const real& val);

    /// Construct a Param object with given value.
    Param(double val);

    /// Construct a Param object with identifier @p id and value @p val.
    Param(const String& id, const real& val);

    /// Return a deep copy of this Param object.
    auto clone() const -> Param;

    /// Assign the attributes of another Param object to this.
    auto assign(const Param& other) -> Param&;

    /// Set the value of the parameter.
    auto value(const real& val) -> Param&;

    /// Return the value of the parameter.
    auto value() const -> const real&;

    /// Return the value of the parameter.
    auto value() -> real&;

    /// Set the unique identifier of the parameter.
    auto id(String id) -> Param&;

    /// Return the unique identifier of the parameter.
    auto id() const -> const String&;

    /// Set the lower bound of the parameter.
    auto lowerbound(double val) -> Param&;

    /// Return the lower bound of the parameter.
    auto lowerbound() const -> double;

    /// Set the upper bound of the parameter.
    auto upperbound(double val) -> Param&;

    /// Return the upper bound of the parameter.
    auto upperbound() const -> double;

    /// Set the parameter to constant or non-constant modes.
    auto isconst(bool val) -> Param&;

    /// Return true if the parameter is constant.
    auto isconst() const -> bool;

    /// Assign a value to this Param object.
    auto operator=(double val) -> Param&;

    /// Assign a value to this Param object.
    auto operator=(const real& val) -> Param&;

    /// Convert this Param object into its value type.
    operator const real&() const;

    /// Convert this Param object into its value type.
    operator real&();

    /// Convert this Param object into its value type.
    operator double() const;

    /// Return a Param object that represents a constant parameter.
    static auto Constant(const real& val) -> Param;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro

//======================================================================
// CODE BELOW NEEDED FOR ARITHMETIC OPERATIONS INVOLVING PARAM
//======================================================================

namespace Reaktoro {

inline auto operator+(const Param& p) { return  p.value(); }
inline auto operator-(const Param& p) { return -p.value(); }

inline auto operator+(const Param& p, const Param& q) { return p.value() + q.value(); }
inline auto operator-(const Param& p, const Param& q) { return p.value() - q.value(); }
inline auto operator*(const Param& p, const Param& q) { return p.value() * q.value(); }
inline auto operator/(const Param& p, const Param& q) { return p.value() / q.value(); }

template<typename T, Requires<isNumeric<T>> = true> auto operator+(const Param& p, const T& x) { return p.value() + x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator-(const Param& p, const T& x) { return p.value() - x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator*(const Param& p, const T& x) { return p.value() * x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator/(const Param& p, const T& x) { return p.value() / x; }

template<typename T, Requires<isNumeric<T>> = true> auto operator+(const T& x, const Param& p) { return x + p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator-(const T& x, const Param& p) { return x - p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator*(const T& x, const Param& p) { return x * p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator/(const T& x, const Param& p) { return x / p.value(); }

} // namespace Reaktoro

//======================================================================
// CODE BELOW NEEDED FOR COMPARISON OPERATIONS INVOLVING PARAM
//======================================================================

namespace Reaktoro {

inline auto operator==(const Param& p, const Param& q) { return p.value() == q.value(); }
inline auto operator!=(const Param& p, const Param& q) { return p.value() != q.value(); }
inline auto operator< (const Param& p, const Param& q) { return p.value()  < q.value(); }
inline auto operator> (const Param& p, const Param& q) { return p.value()  > q.value(); }
inline auto operator<=(const Param& p, const Param& q) { return p.value() <= q.value(); }
inline auto operator>=(const Param& p, const Param& q) { return p.value() >= q.value(); }

template<typename T, Requires<isNumeric<T>> = true> auto operator==(const Param& p, const T& x) { return p.value() == x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator!=(const Param& p, const T& x) { return p.value() != x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator< (const Param& p, const T& x) { return p.value()  < x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator> (const Param& p, const T& x) { return p.value()  > x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator<=(const Param& p, const T& x) { return p.value() <= x; }
template<typename T, Requires<isNumeric<T>> = true> auto operator>=(const Param& p, const T& x) { return p.value() >= x; }

template<typename T, Requires<isNumeric<T>> = true> auto operator==(const T& x, const Param& p) { return x == p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator!=(const T& x, const Param& p) { return x != p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator< (const T& x, const Param& p) { return x  < p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator> (const T& x, const Param& p) { return x  > p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator<=(const T& x, const Param& p) { return x <= p.value(); }
template<typename T, Requires<isNumeric<T>> = true> auto operator>=(const T& x, const Param& p) { return x >= p.value(); }

} // namespace Reaktoro

//======================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING PARAM
//======================================================================

namespace Reaktoro {

template<typename T>
struct MemoizationTraits;

/// Specialize MemoizationTraits for Param.
template<>
struct MemoizationTraits<Param>
{
    using CacheType = real;
};

} // namespace Reaktoro

//======================================================================
// CODE BELOW NEEDED FOR AUTOMATIC DIFFERENTIATION INVOLVING PARAM
//======================================================================

namespace Reaktoro {

template<size_t order, typename U>
auto seed(Param& param, U&& seedval)
{
    autodiff::detail::seed<order>(param.value(), seedval);
}

} // namespace Reaktoro

namespace autodiff {
namespace detail {

/// Implementation of NumberTraits for Reaktoro::Param.
template<>
struct NumberTraits<Reaktoro::Param>
{
    /// The underlying floating point type of Param.
    using NumericType = double;

    /// The order of Param.
    static constexpr auto Order = 1;
};

} // namespace autodiff
} // namespace detail
