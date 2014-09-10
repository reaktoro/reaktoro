// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "ThermoScalar.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>

namespace Reaktor {

ThermoScalar::ThermoScalar()
: m_val(INFINITY), m_ddt(INFINITY), m_ddp(INFINITY), m_ddn()
{}

auto ThermoScalar::val(double val) -> ThermoScalar&
{
    m_val = val;
    return *this;
}

auto ThermoScalar::ddt(double ddt) -> ThermoScalar&
{
    m_ddt = ddt;
    return *this;
}

auto ThermoScalar::ddp(double ddp) -> ThermoScalar&
{
    m_ddp = ddp;
    return *this;
}

auto ThermoScalar::ddn(const Vector& ddn) -> ThermoScalar&
{
    m_ddn = ddn;
    return *this;
}

auto ThermoScalar::val() const -> const double&
{
    if(m_val != INFINITY)
        return m_val;
    else error("Cannot return the value of the thermodynamic quantity.", "Data not initialised.")
}

auto ThermoScalar::ddt() const -> const double&
{
    if(m_ddt != INFINITY)
        return m_ddt;
    else error("Cannot return the partial temperature derivative of the thermodynamic quantity.", "Data not initialised.")
}

auto ThermoScalar::ddp() const -> const double&
{
    if(m_ddp != INFINITY)
        return m_ddp;
    else error("Cannot return the partial pressure derivative of the thermodynamic quantity.", "Data not initialised.")
}

auto ThermoScalar::ddn() const -> const Vector&
{
    if(not m_ddn.empty())
        return m_ddn;
    else error("Cannot return the partial molar derivatives of the thermodynamic quantity.", "Data not initialised.")
}

} // namespace Reaktor
