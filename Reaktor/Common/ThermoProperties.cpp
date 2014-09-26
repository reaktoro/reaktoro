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

#include "ThermoProperties.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/ThermoProperty.hpp>

namespace Reaktor {
namespace internal {

auto assertThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp) -> void
{
    Assert(val.size() == ddt.size() and val.size() == ddp.size(),
        "ThermoProperties requires arguments with the same dimensions.");
}

} // namespace internal

ThermoProperties::ThermoProperties(unsigned nrows)
: m_val(nrows), m_ddt(nrows), m_ddp(nrows)
{}

ThermoProperties::ThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp)
: m_val(val), m_ddt(ddt), m_ddp(ddp)
{
    internal::assertThermoProperties(val, ddt, ddp);
}

auto ThermoProperties::val() const -> const Vector&
{
	return m_val;	
}

auto ThermoProperties::ddt() const -> const Vector&
{
	return m_ddt;	
}

auto ThermoProperties::ddp() const -> const Vector&
{
	return m_ddp;	
}

auto ThermoProperties::row(unsigned irow) -> ThermoPropertiesRow
{
	return ThermoPropertiesRow(*this, irow);
}

auto ThermoProperties::row(unsigned irow) const -> ThermoPropertiesConstRow
{
	return ThermoPropertiesConstRow(*this, irow);
}

ThermoPropertiesRow::ThermoPropertiesRow(const ThermoProperties& properties, unsigned irow)
: val(properties.val().row(irow)),
  ddt(properties.ddt().row(irow)),
  ddp(properties.ddp().row(irow))
{}

ThermoPropertiesConstRow::ThermoPropertiesConstRow(const ThermoProperties& properties, unsigned irow)
: val(properties.val().row(irow)),
  ddt(properties.ddt().row(irow)),
  ddp(properties.ddp().row(irow))
{}

auto ThermoPropertiesRow::operator=(const ThermoProperty& property) -> ThermoPropertiesRow&
{
	val = property.val();
	ddt = property.ddt();
	ddp = property.ddp();
	return *this;
}

auto operator==(const ThermoProperties& l, const ThermoProperties& r) -> bool
{
    return arma::all(l.val() == r.val()) and
           arma::all(l.ddt() == r.ddt()) and
           arma::all(l.ddp() == r.ddp());
}

} // namespace Reaktor
