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

auto ThermoVector::val(const Vector& val) -> ThermoVector&
{
    m_val = val;
    return *this;
}

auto ThermoVector::ddt(const Vector& ddt) -> ThermoVector&
{
    m_ddt = ddt;
    return *this;
}

auto ThermoVector::ddp(const Vector& ddp) -> ThermoVector&
{
    m_ddp = ddp;
    return *this;
}

auto ThermoVector::ddn(const Matrix& ddn) -> ThermoVector&
{
    m_ddn = ddn;
    return *this;
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

auto ThermoVector::ddn() const -> const Matrix&
{
    return m_ddn;
}

auto ThermoVector::row(const Index& irow) -> ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

auto ThermoVector::row(const Index& irow) const -> const ThermoVectorRow
{
    return ThermoVectorRow(*this, irow);
}

ThermoVectorRow::ThermoVectorRow(const ThermoVector& vector, unsigned irow)
: m_val(vector.val().row(irow)), m_ddt(vector.ddt().row(irow)),
  m_ddp(vector.ddp().row(irow)), m_ddn(vector.ddn().row(irow))
{}

auto ThermoVectorRow::operator=(const ThermoScalar& scalar) -> ThermoVectorRow&
{
    m_val = scalar.val();
    m_ddt = scalar.ddt();
    m_ddp = scalar.ddp();
    m_ddn = scalar.ddn();
    return *this;
}

auto ThermoVectorRow::val() const -> const VectorRow&
{
    return m_val;
}

auto ThermoVectorRow::ddt() const -> const VectorRow&
{
    return m_ddt;
}

auto ThermoVectorRow::ddp() const -> const VectorRow&
{
    return m_ddp;
}

auto ThermoVectorRow::ddn() const -> const MatrixRow&
{
    return m_ddn;
}

} // namespace Reaktor
