// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

ThermoVector::ThermoVector(unsigned nrows)
: val(nrows), ddt(nrows), ddp(nrows)
{}

ThermoVector::ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp)
: val(val), ddt(ddt), ddp(ddp)
{
    Assert(val.size() == ddt.size() && val.size() == ddp.size(),
        "Could not construct a ThermoVector instance.",
        "ThermoVector requires arguments with the same dimensions.");
}

auto ThermoVector::row(unsigned irow) -> ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

auto ThermoVector::row(unsigned irow) const -> ThermoVectorConstRow
{
    return ThermoVectorConstRow(*this, irow);
}

auto ThermoVector::rows(unsigned irow, unsigned nrows) -> ThermoVectorRows
{
    return ThermoVectorRows(*this, irow, nrows);
}

auto ThermoVector::rows(unsigned irow, unsigned nrows) const -> ThermoVectorConstRows
{
    return ThermoVectorConstRows(*this, irow, nrows);
}

auto ThermoVector::operator[](unsigned irow) -> ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

auto ThermoVector::operator[](unsigned irow) const -> ThermoVectorConstRow
{
    return ThermoVectorConstRow(*this, irow);
}

ThermoVectorRow::ThermoVectorRow(ThermoVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow])
{}

ThermoVectorConstRow::ThermoVectorConstRow(const ThermoVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow])
{}

ThermoVectorRows::ThermoVectorRows(ThermoVector& vector, unsigned irow, unsigned nrows)
: val(vector.val.segment(irow, nrows)),
  ddt(vector.ddt.segment(irow, nrows)),
  ddp(vector.ddp.segment(irow, nrows))
{}

ThermoVectorConstRows::ThermoVectorConstRows(const ThermoVector& vector, unsigned irow, unsigned nrows)
: val(vector.val.segment(irow, nrows)),
  ddt(vector.ddt.segment(irow, nrows)),
  ddp(vector.ddp.segment(irow, nrows))
{}

auto ThermoVectorRow::operator=(const ThermoScalar& scalar) -> ThermoVectorRow&
{
    val = scalar.val;
    ddt = scalar.ddt;
    ddp = scalar.ddp;
    return *this;
}

auto ThermoVectorRows::operator=(const ThermoVectorRows& block) -> ThermoVectorRows&
{
    val = block.val;
    ddt = block.ddt;
    ddp = block.ddp;
    return *this;
}

auto ThermoVectorRows::operator=(const ThermoVector& vector) -> ThermoVectorRows&
{
    val = vector.val;
    ddt = vector.ddt;
    ddp = vector.ddp;
    return *this;
}

auto operator==(const ThermoVector& l, const ThermoVector& r) -> bool
{
    return l.val == r.val &&
           l.ddt == r.ddt &&
           l.ddp == r.ddp;
}

} // namespace Reaktoro
