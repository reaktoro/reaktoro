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

#include "ChemicalScalar.hpp"

// Reaktor includes
#include <Reaktor/Common/ChemicalVector.hpp>

namespace Reaktor {

ChemicalScalar::ChemicalScalar()
{}

ChemicalScalar::ChemicalScalar(unsigned num_species)
: ChemicalScalar(0.0, 0.0, 0.0, zeros(num_species))
{}

ChemicalScalar::ChemicalScalar(double val, double ddt, double ddp, const Vector& ddn)
: val(val), ddt(ddt), ddp(ddp), ddn(ddn)
{}

ChemicalScalar::ChemicalScalar(const ChemicalVectorRow& row)
{
    *this = row;
}

ChemicalScalar::ChemicalScalar(const ChemicalVectorConstRow& row)
{
    *this = row;
}

auto ChemicalScalar::operator=(const ChemicalVectorRow& row) -> ChemicalScalar&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    ddn = row.ddn;
    return *this;
}

auto ChemicalScalar::operator=(const ChemicalVectorConstRow& row) -> ChemicalScalar&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    ddn = row.ddn;
    return *this;
}

auto operator==(const ChemicalScalar& l, const ChemicalScalar& r) -> bool
{
    return l.val == r.val and
           l.ddt == r.ddt and
           l.ddp == r.ddp and
           l.ddn == r.ddn;
}

} // namespace Reaktor
