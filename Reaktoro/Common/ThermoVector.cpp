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

ThermoVector::ThermoVector()
{}

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

auto ThermoVector::resize(unsigned nrows) -> void
{
    val.resize(nrows);
    ddt.resize(nrows);
    ddp.resize(nrows);
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

auto ThermoVector::operator+=(const ThermoVector& other) -> ThermoVector&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    return *this;
}

auto ThermoVector::operator-=(const ThermoVector& other) -> ThermoVector&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    return *this;
}

auto ThermoVector::operator*=(double scalar) -> ThermoVector&
{
    val *= scalar;
    ddt *= scalar;
    ddp *= scalar;
    return *this;
}

auto ThermoVector::operator/=(double scalar) -> ThermoVector&
{
    *this *= 1.0/scalar;
    return *this;
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

ThermoVectorRow::operator ThermoScalar()
{
    return ThermoScalar(val, ddt, ddp);
}

ThermoVectorConstRow::ThermoVectorConstRow(const ThermoVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow])
{}

ThermoVectorConstRow::operator ThermoScalar()
{
    return ThermoScalar(val, ddt, ddp);
}

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

auto operator+(const ThermoVector& l) -> ThermoVector
{
    return l;
}

auto operator-(const ThermoVector& l) -> ThermoVector
{
    return -1.0 * l;
}

auto operator+(const ThermoVector& l, const ThermoVector& r) -> ThermoVector
{
    ThermoVector res = l;
    res += r;
    return res;
}

auto operator-(const ThermoVector& l, const ThermoVector& r) -> ThermoVector
{
    ThermoVector res = l;
    res -= r;
    return res;
}

auto operator*(double scalar, const ThermoVector& r) -> ThermoVector
{
    ThermoVector res = r;
    res *= scalar;
    return res;
}

auto operator*(const ThermoVector& l, double scalar) -> ThermoVector
{
    return scalar * l;
}

auto operator*(const ThermoScalar& scalar, const ThermoVector& r) -> ThermoVector
{
    ThermoVector res = scalar.val * r;
    res.ddt += scalar.ddt * r.val;
    res.ddp += scalar.ddp * r.val;
    return res;
}

auto operator*(const ThermoVector& l, const ThermoScalar& scalar) -> ThermoVector
{
    return scalar * l;
}

auto operator/(double scalar, const ThermoVector& r) -> ThermoVector
{
    const Vector factor = -scalar/(r.val % r.val);
    ThermoVector res;
    res.val = scalar/r.val;
    res.ddt = factor % r.ddt;
    res.ddp = factor % r.ddp;
    return res;
}

auto operator/(const ThermoVector& l, double scalar) -> ThermoVector
{
    return (1.0/scalar) * l;
}

auto operator/(const ThermoScalar& l, const ThermoVector& r) -> ThermoVector
{
    const Vector factor = 1.0/(r.val % r.val);
    ThermoVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) % factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) % factor;
    return res;
}

auto operator/(const ThermoVector& l, const ThermoScalar& r) -> ThermoVector
{
    const double factor = 1.0/(r.val * r.val);
    ThermoVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) * factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) * factor;
    return res;
}

auto operator/(const ThermoVector& l, const ThermoVector& r) -> ThermoVector
{
    const Vector factor = 1.0/(r.val % r.val);
    ThermoVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt % r.val - r.ddt % l.val) % factor;
    res.ddp = (l.ddp % r.val - r.ddp % l.val) % factor;
    return res;
}

auto operator%(const ThermoVector& l, const ThermoVector& r) -> ThermoVector
{
    ThermoVector res;
    res.val = diag(l.val) * r.val;
    res.ddt = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddp = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    return res;
}

auto pow(const ThermoVector& a, double power) -> ThermoVector
{
    ThermoVector b;
    b.val = pow(a.val, power);
    b.ddt = power * diag(b.val) * a.ddt/a.val;
    b.ddp = power * diag(b.val) * a.ddp/a.val;
    return b;
}

auto exp(const ThermoVector& a) -> ThermoVector
{
    ThermoVector b;
    b.val = exp(a.val);
    b.ddt = b.val % a.ddt;
    b.ddp = b.val % a.ddp;
    return b;
}

auto log(const ThermoVector& a) -> ThermoVector
{
    ThermoVector b;
    b.val = log(a.val);
    b.ddt = a.ddt/a.val;
    b.ddp = a.ddp/a.val;
    return b;
}

auto log10(const ThermoVector& a) -> ThermoVector
{
    const double ln10 = 2.302585092994046;
    ThermoVector b = log(a);
    b /= ln10;
    return b;
}

auto sum(const ThermoVector& l) -> ThermoScalar
{
    ThermoScalar res;
    res.val = l.val.sum();
    res.ddt = l.ddt.sum();
    res.ddp = l.ddp.sum();
    return res;
}

} // namespace Reaktoro
