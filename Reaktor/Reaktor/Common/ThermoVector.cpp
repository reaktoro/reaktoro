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
#include <Reaktor/Common/ThermoScalar.hpp>

namespace Reaktor {

ThermoVector::ThermoVector()
{}

ThermoVector::ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn)
: val(val), ddt(ddt), ddp(ddp), ddn(ddn)
{}

ThermoVector::ThermoVector(unsigned nrows, unsigned nspecies)
: val(nrows), ddt(nrows), ddp(nrows), ddn(nrows, nspecies)
{}

auto ThermoVector::row(const Index& irow) -> ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

auto ThermoVector::row(const Index& irow) const -> const ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

auto ThermoVector::zero(unsigned nrows, unsigned nspecies) -> ThermoVector
{
    ThermoVector vector;
    vector.val = zeros(nrows);
    vector.ddt = zeros(nrows);
    vector.ddp = zeros(nrows);
    vector.ddn = zeros(nrows, nspecies);
    return vector;
}

ThermoVectorRow::ThermoVectorRow(const ThermoVector& vector, unsigned irow)
: val(vector.val.row(irow)), ddt(vector.ddt.row(irow)),
  ddp(vector.ddp.row(irow)), ddn(vector.ddn.row(irow))
{}

auto ThermoVectorRow::operator=(const ThermoScalar& scalar) -> ThermoVectorRow&
{
    val = scalar.val;
    ddt = scalar.ddt;
    ddp = scalar.ddp;
    ddn = scalar.ddn;
    return *this;
}

} // namespace Reaktor
