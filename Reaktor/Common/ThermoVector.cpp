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

#include "ThermoVector.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>

namespace Reaktor {
namespace internal {

auto assertThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp) -> void
{
    Assert(val.size() == ddt.size() and val.size() == ddp.size(),
        "ThermoVector requires arguments with the same dimensions.");
}

} // namespace internal

ThermoVector::ThermoVector(unsigned nrows)
: m_val(nrows), m_ddt(nrows), m_ddp(nrows)
{}

ThermoVector::ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp)
: m_val(val), m_ddt(ddt), m_ddp(ddp)
{
    internal::assertThermoVector(val, ddt, ddp);
}

auto ThermoVector::val() const -> const Vector&
{
	return m_val;	
}

auto ThermoVector::ddt() const -> const Vector&
{
	return m_ddt;	
}

auto ThermoVector::ddp() const -> const Vector&
{
	return m_ddp;	
}

auto ThermoVector::row(unsigned irow) -> ThermoVectorRow
{
	return ThermoVectorRow(*this, irow);
}

auto ThermoVector::row(unsigned irow) const -> ThermoVectorConstRow
{
	return ThermoVectorConstRow(*this, irow);
}

ThermoVectorRow::ThermoVectorRow(const ThermoVector& properties, unsigned irow)
: val(properties.val().row(irow)),
  ddt(properties.ddt().row(irow)),
  ddp(properties.ddp().row(irow))
{}

ThermoVectorConstRow::ThermoVectorConstRow(const ThermoVector& properties, unsigned irow)
: val(properties.val().row(irow)),
  ddt(properties.ddt().row(irow)),
  ddp(properties.ddp().row(irow))
{}

auto ThermoVectorRow::operator=(const ThermoScalar& property) -> ThermoVectorRow&
{
	val = property.val();
	ddt = property.ddt();
	ddp = property.ddp();
	return *this;
}

auto operator==(const ThermoVector& l, const ThermoVector& r) -> bool
{
    return arma::all(l.val() == r.val()) and
           arma::all(l.ddt() == r.ddt()) and
           arma::all(l.ddp() == r.ddp());
}

} // namespace Reaktor
