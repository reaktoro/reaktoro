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

#include "ThermoProperty.hpp"

// Reaktor includes
#include <Reaktor/Common/ThermoProperties.hpp>

namespace Reaktor {

ThermoProperty::ThermoProperty()
{}

ThermoProperty::ThermoProperty(double val, double ddt, double ddp)
: m_val(val), m_ddt(ddt), m_ddp(ddp)
{}

auto ThermoProperty::val() const -> double
{
    return m_val;
}

auto ThermoProperty::ddt() const -> double
{
    return m_ddt;
}

auto ThermoProperty::ddp() const -> double
{
    return m_ddp;
}

auto ThermoProperty::operator=(const ThermoPropertiesRow& row) -> ThermoProperty&
{
    m_val = row.val[0];
    m_ddt = row.ddt[0];
    m_ddp = row.ddp[0];
    return *this;
}

auto ThermoProperty::operator=(const ThermoPropertiesConstRow& row) -> ThermoProperty&
{
    m_val = row.val[0];
    m_ddt = row.ddt[0];
    m_ddp = row.ddp[0];
    return *this;
}

}  // namespace Reaktor
