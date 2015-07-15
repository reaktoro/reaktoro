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

#include "ChemicalVector.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {
namespace {

auto assertDimensions(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn) -> void
{
    Assert(val.size() == ddt.size() &&  val.size() == ddt.size() && val.size() == ddn.rows(),
        "Could not construct a ChemicalVector instance.",
        "ChemicalVector requires arguments with the same row-dimensions.");
}

} // namespace

auto ChemicalVector::Composition(const Vector& n) -> ChemicalVector
{
    const unsigned size = n.rows();
    ChemicalVector res(size);
    res.val = n;
    res.ddn = identity(size, size);
    return res;
}

auto ChemicalVector::Composition(std::initializer_list<double> n) -> ChemicalVector
{
    const unsigned size = n.size();
    ChemicalVector res(size);
    res.ddn = identity(size, size);
    unsigned i = 0;
    for(double ni : n)
        res.val[i++] = ni;
    return res;
}

ChemicalVector::ChemicalVector()
{}

ChemicalVector::ChemicalVector(unsigned nspecies)
: ChemicalVector(nspecies, nspecies)
{}

ChemicalVector::ChemicalVector(unsigned nrows, unsigned nspecies)
: val(zeros(nrows)), ddt(zeros(nrows)), ddp(zeros(nrows)), ddn(zeros(nrows, nspecies))
{}

ChemicalVector::ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn)
: val(val), ddt(ddt), ddp(ddp), ddn(ddn)
{
    assertDimensions(val, ddt, ddp, ddn);
}

ChemicalVector::ChemicalVector(const ChemicalVectorRows& rows)
: val(rows.val), ddt(rows.ddt), ddp(rows.ddp), ddn(rows.ddn)
{
    assertDimensions(val, ddt, ddp, ddn);
}

ChemicalVector::ChemicalVector(const ChemicalVectorRowsConst& rows)
: val(rows.val), ddt(rows.ddt), ddp(rows.ddp), ddn(rows.ddn)
{
    assertDimensions(val, ddt, ddp, ddn);
}

ChemicalVector::ChemicalVector(const ChemicalVectorRowsCols& block)
: val(block.val), ddt(block.ddt), ddp(block.ddp), ddn(block.ddn)
{
    assertDimensions(val, ddt, ddp, ddn);
}

ChemicalVector::ChemicalVector(const ChemicalVectorRowsColsConst& block)
: val(block.val), ddt(block.ddt), ddp(block.ddp), ddn(block.ddn)
{
    assertDimensions(val, ddt, ddp, ddn);
}

auto ChemicalVector::resize(unsigned nrows) -> void
{
    resize(nrows, nrows);
}

auto ChemicalVector::resize(unsigned nrows, unsigned ncols) -> void
{
    val.resize(nrows);
    ddt.resize(nrows);
    ddp.resize(nrows);
    ddn.resize(nrows, ncols);
}

auto ChemicalVector::row(unsigned irow) -> ChemicalVectorRow
{
    return ChemicalVectorRow(*this, irow);
}

auto ChemicalVector::row(unsigned irow, unsigned icol, unsigned ncols) -> ChemicalVectorRow
{
    return ChemicalVectorRow(*this, irow, icol, ncols);
}

auto ChemicalVector::row(unsigned irow) const -> ChemicalVectorRowConst
{
    return ChemicalVectorRowConst(*this, irow);
}

auto ChemicalVector::row(unsigned irow, unsigned icol, unsigned ncols) const -> ChemicalVectorRowConst
{
    return ChemicalVectorRowConst(*this, irow, icol, ncols);
}

auto ChemicalVector::rows(unsigned irow, unsigned nrows) -> ChemicalVectorBlock
{
    return ChemicalVectorBlock(*this, irow, nrows);
}

auto ChemicalVector::rows(unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) -> ChemicalVectorBlock
{
    return ChemicalVectorBlock(*this, irow, icol, nrows, ncols);
}

auto ChemicalVector::rows(unsigned irow, unsigned nrows) const -> ChemicalVectorBlockConst
{
    return ChemicalVectorBlockConst(*this, irow, nrows);
}

auto ChemicalVector::rows(unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) const -> ChemicalVectorBlockConst
{
    return ChemicalVectorBlockConst(*this, irow, icol, nrows, ncols);
}

auto ChemicalVector::rows(const Indices& irows) -> ChemicalVectorRows
{
    return ChemicalVectorRows(*this, irows);
}

auto ChemicalVector::rows(const Indices& irows) const -> ChemicalVectorRowsConst
{
    return ChemicalVectorRowsConst(*this, irows);
}

auto ChemicalVector::rows(const Indices& irows, const Indices& icols) -> ChemicalVectorRowsCols
{
    return ChemicalVectorRowsCols(*this, irows, icols);
}

auto ChemicalVector::rows(const Indices& irows, const Indices& icols) const -> ChemicalVectorRowsColsConst
{
    return ChemicalVectorRowsColsConst(*this, irows, icols);
}

auto ChemicalVector::operator+=(const ChemicalVector& other) -> ChemicalVector&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    ddn += other.ddn;
    return *this;
}

auto ChemicalVector::operator+=(const ThermoVector& other) -> ChemicalVector&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    return *this;
}

auto ChemicalVector::operator-=(const ChemicalVector& other) -> ChemicalVector&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    ddn -= other.ddn;
    return *this;
}

auto ChemicalVector::operator-=(const ThermoVector& other) -> ChemicalVector&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    return *this;
}

auto ChemicalVector::operator*=(double scalar) -> ChemicalVector&
{
    val *= scalar;
    ddt *= scalar;
    ddp *= scalar;
    ddn *= scalar;
    return *this;
}

auto ChemicalVector::operator/=(double scalar) -> ChemicalVector&
{
    *this *= 1.0/scalar;
    return *this;
}

auto ChemicalVector::operator[](unsigned irow) -> ChemicalVectorRow
{
    return ChemicalVectorRow(*this, irow);
}

auto ChemicalVector::operator[](unsigned irow) const -> ChemicalVectorRowConst
{
    return ChemicalVectorRowConst(*this, irow);
}

ChemicalVectorRow::ChemicalVectorRow(ChemicalVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow).segment(0, vector.ddn.cols()))
{}

ChemicalVectorRow::ChemicalVectorRow(ChemicalVector& vector, unsigned irow, unsigned icol, unsigned ncols)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow).segment(icol, ncols))
{}

ChemicalVectorRowConst::ChemicalVectorRowConst(const ChemicalVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow).segment(0, vector.ddn.cols()))
{}

ChemicalVectorRowConst::ChemicalVectorRowConst(const ChemicalVector& vector, unsigned irow, unsigned icol, unsigned ncols)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow).segment(icol, ncols))
{}

ChemicalVectorRows::ChemicalVectorRows(ChemicalVector& vector, const Indices& irows)
: val(rows(vector.val, irows)),
  ddt(rows(vector.ddt, irows)),
  ddp(rows(vector.ddp, irows)),
  ddn(rows(vector.ddn, irows))
{}

ChemicalVectorRowsConst::ChemicalVectorRowsConst(const ChemicalVector& vector, const Indices& irows)
: val(rows(vector.val, irows)),
  ddt(rows(vector.ddt, irows)),
  ddp(rows(vector.ddp, irows)),
  ddn(rows(vector.ddn, irows))
{}

ChemicalVectorRowsCols::ChemicalVectorRowsCols(ChemicalVector& vector, const Indices& irows, const Indices& icols)
: val(rows(vector.val, irows)),
  ddt(rows(vector.ddt, irows)),
  ddp(rows(vector.ddp, irows)),
  ddn(submatrix(vector.ddn, irows, icols))
{}

ChemicalVectorRowsColsConst::ChemicalVectorRowsColsConst(const ChemicalVector& vector, const Indices& irows, const Indices& icols)
: val(rows(vector.val, irows)),
  ddt(rows(vector.ddt, irows)),
  ddp(rows(vector.ddp, irows)),
  ddn(submatrix(vector.ddn, irows, icols))
{}

auto ChemicalVectorRows::operator=(const ChemicalVector& vector) -> ChemicalVectorRows&
{
    val = vector.val;
    ddt = vector.ddt;
    ddp = vector.ddp;
    ddn = vector.ddn;
    return *this;
}

auto ChemicalVectorRowsCols::operator=(const ChemicalVector& vector) -> ChemicalVectorRowsCols&
{
    val = vector.val;
    ddt = vector.ddt;
    ddp = vector.ddp;
    ddn = vector.ddn;
    return *this;
}

ChemicalVectorBlock::ChemicalVectorBlock(ChemicalVector& vector, unsigned irow, unsigned nrows)
: ChemicalVectorBlock(vector, irow, 0, nrows, vector.ddn.cols())
{}

ChemicalVectorBlock::ChemicalVectorBlock(ChemicalVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols)
: val(vector.val.segment(irow, nrows)),
  ddt(vector.ddt.segment(irow, nrows)),
  ddp(vector.ddp.segment(irow, nrows)),
  ddn(vector.ddn.block(irow, icol, nrows, ncols))
{}

ChemicalVectorBlockConst::ChemicalVectorBlockConst(const ChemicalVector& vector, unsigned irow, unsigned nrows)
: ChemicalVectorBlockConst(vector, irow, 0, nrows, vector.ddn.cols())
{}

ChemicalVectorBlockConst::ChemicalVectorBlockConst(const ChemicalVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols)
: val(vector.val.segment(irow, nrows)),
  ddt(vector.ddt.segment(irow, nrows)),
  ddp(vector.ddp.segment(irow, nrows)),
  ddn(vector.ddn.block(irow, icol, nrows, ncols))
{}

auto ChemicalVectorBlock::operator=(const ChemicalVectorBlock& block) -> ChemicalVectorBlock&
{
    val = block.val;
    ddt = block.ddt;
    ddp = block.ddp;
    ddn = block.ddn;
    return *this;
}

auto ChemicalVectorBlock::operator=(const ChemicalVector& vector) -> ChemicalVectorBlock&
{
    val = vector.val;
    ddt = vector.ddt;
    ddp = vector.ddp;
    ddn = vector.ddn;
    return *this;
}

auto ChemicalVectorRow::operator=(const ChemicalVectorRow& row) -> ChemicalVectorRow&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    ddn = row.ddn;
    return *this;
}

auto ChemicalVectorRow::operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&
{
    val = scalar.val;
    ddt = scalar.ddt;
    ddp = scalar.ddp;
    ddn = scalar.ddn;
    return *this;
}

auto ChemicalVectorRow::operator+=(const ChemicalVectorRow& other) -> ChemicalVectorRow&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    ddn += other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator+=(const ChemicalVectorRowConst& other) -> ChemicalVectorRow&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    ddn += other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator+=(const ChemicalScalar& other) -> ChemicalVectorRow&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    ddn += other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator+=(const ThermoScalar& other) -> ChemicalVectorRow&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    return *this;
}

auto ChemicalVectorRow::operator-=(const ChemicalVectorRow& other) -> ChemicalVectorRow&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    ddn -= other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator-=(const ChemicalVectorRowConst& other) -> ChemicalVectorRow&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    ddn -= other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator-=(const ChemicalScalar& other) -> ChemicalVectorRow&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    ddn -= other.ddn;
    return *this;
}

auto ChemicalVectorRow::operator-=(const ThermoScalar& other) -> ChemicalVectorRow&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    return *this;
}


auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool
{
    return l.val == r.val &&
           l.ddt == r.ddt &&
           l.ddp == r.ddp &&
           l.ddn == r.ddn;
}

auto operator+(const ChemicalVector& l) -> ChemicalVector
{
    return l;
}

auto operator-(const ChemicalVector& l) -> ChemicalVector
{
    return -1.0 * l;
}

auto operator+(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = l;
    res += r;
    return res;
}

auto operator+(const ChemicalVector& l, const ThermoVector& r) -> ChemicalVector
{
    ChemicalVector res = l;
    res += r;
    return res;
}

auto operator+(const ThermoVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = r;
    res += l;
    return res;
}

auto operator-(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = l;
    res -= r;
    return res;
}

auto operator-(const ChemicalVector& l, const ThermoVector& r) -> ChemicalVector
{
    ChemicalVector res = l;
    res -= r;
    return res;
}

auto operator-(const ThermoVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = r;
    res -= l;
    return res;
}

auto operator*(double scalar, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = r;
    res *= scalar;
    return res;
}

auto operator*(const ChemicalVector& l, double scalar) -> ChemicalVector
{
    return scalar * l;
}

auto operator*(const ThermoScalar& scalar, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res = scalar.val * r;
    res.ddt += scalar.ddt * r.val;
    res.ddp += scalar.ddp * r.val;
    return res;
}

auto operator*(const ChemicalVector& l, const ThermoScalar& scalar) -> ChemicalVector
{
    return scalar * l;
}

auto operator/(double scalar, const ChemicalVector& r) -> ChemicalVector
{
    const Vector factor = -scalar/(r.val % r.val);
    ChemicalVector res;
    res.val = scalar/r.val;
    res.ddt = factor % r.ddt;
    res.ddp = factor % r.ddp;
    res.ddn = diag(factor) * r.ddn;
    return res;
}

auto operator/(const ChemicalVector& l, double scalar) -> ChemicalVector
{
    return (1.0/scalar) * l;
}

auto operator/(const ThermoScalar& l, const ChemicalVector& r) -> ChemicalVector
{
    const Vector factor = 1.0/(r.val % r.val);
    ChemicalVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) % factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) % factor;
    res.ddn = diag(factor) * (- l.val * r.ddn);
    return res;
}

auto operator/(const ChemicalVector& l, const ThermoScalar& r) -> ChemicalVector
{
    const double factor = 1.0/(r.val * r.val);
    ChemicalVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) * factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) * factor;
    res.ddn = factor * (r.val * l.ddn);
    return res;
}

auto operator/(const ChemicalScalar& l, const ChemicalVector& r) -> ChemicalVector
{
    const Vector factor = 1.0/(r.val % r.val);
    ChemicalVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) % factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) % factor;
    res.ddn = diag(factor) * (r.val * tr(l.ddn) - l.val * r.ddn);
    return res;
}

auto operator/(const ChemicalVector& l, const ChemicalScalar& r) -> ChemicalVector
{
    const double factor = 1.0/(r.val * r.val);
    ChemicalVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt * r.val - r.ddt * l.val) * factor;
    res.ddp = (l.ddp * r.val - r.ddp * l.val) * factor;
    res.ddn = factor * (r.val * l.ddn - l.val * tr(r.ddn));
    return res;
}

auto operator/(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector
{
    const Vector factor = 1.0/(r.val % r.val);
    ChemicalVector res;
    res.val = l.val / r.val;
    res.ddt = (l.ddt % r.val - r.ddt % l.val) % factor;
    res.ddp = (l.ddp % r.val - r.ddp % l.val) % factor;
    res.ddn = diag(factor) * (diag(r.val) * l.ddn - diag(l.val) * r.ddn);
    return res;
}

auto operator%(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res;
    res.val = diag(l.val) * r.val;
    res.ddt = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddp = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddn = diag(l.val) * r.ddn + diag(r.val) * l.ddn;
    return res;
}

auto operator%(const ThermoVector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res;
    res.val = diag(l.val) * r.val;
    res.ddt = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddp = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddn = diag(l.val) * r.ddn;
    return res;
}

auto operator%(const ChemicalVector& l, const ThermoVector& r) -> ChemicalVector
{
    ChemicalVector res;
    res.val = diag(l.val) * r.val;
    res.ddt = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddp = diag(l.val) * r.ddt + diag(r.val) * l.ddt;
    res.ddn = diag(r.val) * l.ddn;
    return res;
}

auto operator%(const Vector& l, const ChemicalVector& r) -> ChemicalVector
{
    ChemicalVector res;
    res.val = diag(l) * r.val;
    res.ddt = diag(l) * r.ddt;
    res.ddp = diag(l) * r.ddp;
    res.ddn = diag(l) * r.ddn;
    return res;
}

auto operator%(const ChemicalVector& l, const Vector& r) -> ChemicalVector
{
    return r % l;
}

auto pow(const ChemicalVector& a, double power) -> ChemicalVector
{
    ChemicalVector b;
    b.val = pow(a.val, power);
    b.ddt = power * diag(b.val) * a.ddt/a.val;
    b.ddp = power * diag(b.val) * a.ddp/a.val;
    b.ddn = power * diag(b.val) * a.ddn/a.val;
    return b;
}

auto exp(const ChemicalVector& a) -> ChemicalVector
{
    ChemicalVector b;
    b.val = exp(a.val);
    b.ddt = b.val % a.ddt;
    b.ddp = b.val % a.ddp;
    b.ddn = diag(b.val) * a.ddn;
    return b;
}

auto log(const ChemicalVector& a) -> ChemicalVector
{
    ChemicalVector b;
    b.val = log(a.val);
    b.ddt = a.ddt/a.val;
    b.ddp = a.ddp/a.val;
    b.ddn = diag(inv(a.val)) * a.ddn;
    return b;
}

auto log10(const ChemicalVector& a) -> ChemicalVector
{
    const double ln10 = 2.302585092994046;
    ChemicalVector b = log(a);
    b /= ln10;
    return b;
}

auto sum(const ChemicalVector& l) -> ChemicalScalar
{
    ChemicalScalar res;
    res.val = l.val.sum();
    res.ddt = l.ddt.sum();
    res.ddp = l.ddp.sum();
    res.ddn = l.ddn.rowwise().sum();
    return res;
}

} // namespace Reaktoro
