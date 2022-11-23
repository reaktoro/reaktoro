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

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/NumberTraits.hpp>

namespace Reaktoro {

/// The result of a reaction rate model evaluation.
class ReactionRate
{
public:
    /// Construct a default ReactionRate object.
    ReactionRate()
    {}

    /// Construct a ReactionRate object with given rate value.
    template<typename T, Requires<isNumeric<T>> = true>
    ReactionRate(T const& value)
    : m_value(value) {}

    /// Return a ReactionRate object that represents the residual of an enforced equation `f(props) = 0` instead of a reaction rate.
    static auto enforce(real const& value) -> ReactionRate
    {
        ReactionRate res(value);
        res.m_equation_mode = true;
        return res;
    }

    /// Return the underlying real object in the ReactionRate object.
    auto value() const -> real const&
    {
        return m_value;
    }

    /// Return true if this ReactionRate object is in equation mode enabled by @ref enforce.
    auto onEquationMode() const -> bool
    {
        return m_equation_mode;
    }

    /// Convert this ReactionRate object into a real object.
    operator real const&() const
    {
        return m_value;
    }

    /// Assign a value to this ReactionRate object.
    template<typename T, Requires<isNumeric<T>> = true>
    auto operator=(T const& value) -> ReactionRate&
    {
        m_value = value;
        return *this;
    }

    template<typename T, Requires<isNumeric<T>> = true> auto operator+=(T const& scalar) -> ReactionRate& { m_value += scalar; return *this; }
    template<typename T, Requires<isNumeric<T>> = true> auto operator-=(T const& scalar) -> ReactionRate& { m_value -= scalar; return *this; }
    template<typename T, Requires<isNumeric<T>> = true> auto operator*=(T const& scalar) -> ReactionRate& { m_value *= scalar; return *this; }
    template<typename T, Requires<isNumeric<T>> = true> auto operator/=(T const& scalar) -> ReactionRate& { m_value /= scalar; return *this; }

private:
    /// The computed value of the reaction rate or the residual of an equation that the rate model is enforcing.
    real m_value = {};

    /// The boolean flag that indices whether `value` is to be interpreted as the residual of an enforced equation `f(props) = 0` instead of a reaction rate.
    bool m_equation_mode = false;
};

inline auto operator+(ReactionRate const& rate) { return rate; }
inline auto operator-(ReactionRate rate) { return rate *= -1.0; }

template<typename T, Requires<isNumeric<T>> = true> auto operator+(ReactionRate rate, T const& scalar) { return rate += scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator-(ReactionRate rate, T const& scalar) { return rate -= scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator*(ReactionRate rate, T const& scalar) { return rate *= scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator/(ReactionRate rate, T const& scalar) { return rate /= scalar; }

template<typename T, Requires<isNumeric<T>> = true> auto operator+(T const& scalar, ReactionRate rate) { return rate += scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator-(T const& scalar, ReactionRate rate) { return rate -= scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator*(T const& scalar, ReactionRate rate) { return rate *= scalar; }
template<typename T, Requires<isNumeric<T>> = true> auto operator/(T const& scalar, ReactionRate rate) { return rate /= scalar; }

} // namespace Reaktoro
