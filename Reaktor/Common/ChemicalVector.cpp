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
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>

namespace Reaktor {

ChemicalVector::ChemicalVector()
{}

ChemicalVector::ChemicalVector(unsigned nrows, unsigned ncols)
: val(zeros(nrows)), ddt(zeros(nrows)), ddp(zeros(nrows)), ddn(zeros(nrows, ncols))
{}

ChemicalVector::ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn)
: val(val), ddt(ddt), ddp(ddp), ddn(ddn)
{
    Assert(val.size() == ddt.size() and val.size() == ddt.size() and val.size() == ddn.rows(),
        "Could not construct a ChemicalVector instance.",
        "ChemicalVector requires arguments with the same dimensions.");
}

auto ChemicalVector::row(unsigned irow) -> ChemicalVectorRow
{
    return ChemicalVectorRow(*this, irow);
}

auto ChemicalVector::row(unsigned irow) const -> ChemicalVectorConstRow
{
    return ChemicalVectorConstRow(*this, irow);
}

ChemicalVectorRow::ChemicalVectorRow(ChemicalVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow))
{}

ChemicalVectorConstRow::ChemicalVectorConstRow(const ChemicalVector& vector, unsigned irow)
: val(vector.val[irow]),
  ddt(vector.ddt[irow]),
  ddp(vector.ddp[irow]),
  ddn(vector.ddn.row(irow))
{}

auto ChemicalVectorRow::operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&
{
    val = scalar.val();
    ddt = scalar.ddt();
    ddp = scalar.ddp();
    ddn = scalar.ddn();
    return *this;
}

auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool
{
    return l.val() == r.val() and
           l.ddt() == r.ddt() and
           l.ddp() == r.ddp() and
           l.ddn() == r.ddn();
}

} // namespace Reaktor
