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

namespace Reaktoro {

/// The result of a reaction rate model evaluation.
class Rate
{
public:
    /// Construct a default Rate object.
    Rate()
    {}

    /// Construct a Rate object with given rate value.
    Rate(real const& value)
    : m_value(m_value) {}

    /// Return a Rate object that represents the residual of an enforced equation `f(props) = 0` instead of a reaction rate.
    static Rate enforce(real const& value)
    {
        Rate res(value);
        res.m_equation_mode = true;
        return res;
    }

    /// Convert this Rate object into a real object.
    operator real const&() const
    {
        return m_value;
    }

    /// Assign a real value to this Rate object.
    auto operator=(real const& value) -> Rate&
    {
        m_value = value;
        return *this;
    }

private:
    /// The computed value of the reaction rate or the residual of an equation that the rate model is enforcing.
    real m_value = {};

    /// The boolean flag that indices whether `value` is to be interpreted as the residual of an enforced equation `f(props) = 0` instead of a reaction rate.
    bool m_equation_mode = false;
};

} // namespace Reaktoro
