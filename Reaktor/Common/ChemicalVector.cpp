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

#include "ChemicalVector.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>

namespace Reaktor {
namespace {

auto assertChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn) -> void
{
    Assert(val.size() == ddt.size() and val.size() == ddt.size() and val.size() == ddn.n_rows,
        "ChemicalVector requires arguments with the same dimensions.");
}

} // namespace

ChemicalVector::ChemicalVector()
{}

ChemicalVector::ChemicalVector(unsigned nrows, unsigned ncols)
: m_val(nrows), m_ddt(nrows), m_ddp(nrows), m_ddn(nrows, ncols)
{}

ChemicalVector::ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn)
: m_val(val), m_ddt(ddt), m_ddp(ddp), m_ddn(ddn)
{
    assertChemicalVector(val, ddt, ddp, ddn);
}

auto ChemicalVector::val() const -> const Vector&
{
    return m_val;
}

auto ChemicalVector::ddt() const -> const Vector&
{
    return m_ddt;
}

auto ChemicalVector::ddp() const -> const Vector&
{
    return m_ddp;
}

auto ChemicalVector::ddn() const -> const Matrix&
{
    return m_ddn;
}

auto ChemicalVector::row(unsigned irow) -> ChemicalVectorRow
{
	return ChemicalVectorRow(*this, irow);
}

auto ChemicalVector::row(unsigned irow) const -> ChemicalVectorConstRow
{
	return ChemicalVectorConstRow(*this, irow);
}

ChemicalVectorRow::ChemicalVectorRow(const ChemicalVector& vector, unsigned irow)
: val(vector.val().row(irow)),
  ddt(vector.ddt().row(irow)),
  ddp(vector.ddp().row(irow)),
  ddn(vector.ddn().row(irow))
{}

ChemicalVectorConstRow::ChemicalVectorConstRow(const ChemicalVector& vector, unsigned irow)
: val(vector.val().row(irow)),
  ddt(vector.ddt().row(irow)),
  ddp(vector.ddp().row(irow)),
  ddn(vector.ddn().row(irow))
{}

auto ChemicalVectorRow::operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&
{
	val = scalar.val();
	ddt = scalar.ddt();
	ddp = scalar.ddp();
	ddn = scalar.ddn().t();
	return *this;
}

auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool
{
    return arma::all(l.val() == r.val()) and
           arma::all(l.ddt() == r.ddt()) and
           arma::all(l.ddp() == r.ddp()) and
           arma::all(arma::all(l.ddn() == r.ddn()));
}

} // namespace Reaktor
