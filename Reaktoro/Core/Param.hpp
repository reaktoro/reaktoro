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
#include <limits>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>

namespace Reaktoro {

/// A type used to represent the value of a parameter and its lower and upper bounds.
class Param
{
public:
    /// Construct a default Param object.
    Param() {}

    /// Construct a Param object with given value.
    Param(const real& val) : m_value(val) {}

    /// Construct a Param object with given value.
    template<typename T, EnableIf<isArithmetic<T>>...>
    Param(const T& val) : m_value(val) {}

    /// Set the value of the parameter.
    auto value(const real& val) -> Param& { warning(val <= m_lowerbound || val >= m_upperbound, "Setting parameter with value ", val, " violates either its lower bound (", m_lowerbound, ") or its upper bound (", m_upperbound, ")."); m_value = val; return *this; }

    /// Return the value of the parameter.
    auto value() const -> const real& { return m_value; }

    /// Set the lower bound of the parameter.
    auto lowerbound(double val) -> Param& { m_lowerbound = val; return *this; }

    /// Return the lower bound of the parameter.
    auto lowerbound() const -> double { return m_lowerbound; }

    /// Set the upper bound of the parameter.
    auto upperbound(double val) -> Param& { m_upperbound = val; return *this; }

    /// Return the upper bound of the parameter.
    auto upperbound() const -> double { return m_upperbound; }

    /// Set the parameter to constant or non-constant modes.
    auto isconst(bool val) -> Param& { m_isconst = val; return *this; }

    /// Return true if the parameter is constant.
    auto isconst() const -> bool { return m_isconst; }

    /// Assign a real value to this parameter value.
    auto operator=(const real& val) -> Param& { value(val); return *this; }

    /// Convert this Param object into its value type.
    operator const real&() const { return m_value; }

    /// Return a Param object that represents a constant parameter.
    static auto Constant(const real& val) -> Param { Param param(val); param.isconst(true); return param; }

private:
    /// The parameter value.
    real m_value = {};

    /// The parameter lower bound value (default: -inf).
    double m_lowerbound = -std::numeric_limits<double>::infinity();

    /// The parameter upper bound value (default: +inf).
    double m_upperbound = +std::numeric_limits<double>::infinity();

    /// The boolean flag that indicates if this parameter is constant.
    double m_isconst = false;
};

} // namespace Reaktoro
