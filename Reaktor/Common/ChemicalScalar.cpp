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

ChemicalScalar::ChemicalScalar(double val, double ddt, double ddp, const Vector& ddn)
: m_val(val), m_ddt(ddt), m_ddp(ddp), m_ddn(ddn)
{}

ChemicalScalar::ChemicalScalar(const ChemicalVectorRow& row)
{
    *this = row;
}

ChemicalScalar::ChemicalScalar(const ChemicalVectorConstRow& row)
{
    *this = row;
}

auto ChemicalScalar::val() const -> double
{
    return m_val;
}

auto ChemicalScalar::ddt() const -> double
{
    return m_ddt;
}

auto ChemicalScalar::ddp() const -> double
{
    return m_ddp;
}

auto ChemicalScalar::ddn() const -> const Vector&
{
    return m_ddn;
}

auto ChemicalScalar::operator=(const ChemicalVectorRow& row) -> ChemicalScalar&
{
    m_val = row.val[0];
    m_ddt = row.ddt[0];
    m_ddp = row.ddp[0];
    m_ddn = row.ddn.t();
    return *this;
}

auto ChemicalScalar::operator=(const ChemicalVectorConstRow& row) -> ChemicalScalar&
{
    m_val = row.val[0];
    m_ddt = row.ddt[0];
    m_ddp = row.ddp[0];
    m_ddn = row.ddn.t();
    return *this;
}

auto operator==(const ChemicalScalar& l, const ChemicalScalar& r) -> bool
{
    return l.val() == r.val() and
           l.ddt() == r.ddt() and
           l.ddp() == r.ddp() and
           arma::all(l.ddn() == r.ddn());
}

} // namespace Reaktor
