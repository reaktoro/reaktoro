// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "ChemicalField.hpp"

namespace Reaktoro {

ChemicalField::ChemicalField()
{}

auto ChemicalField::resize(unsigned num_points, unsigned num_components) -> void
{
    m_val.resize(num_points);

    if(m_ddt.size()) m_ddt.resize(num_points);
    if(m_ddp.size()) m_ddp.resize(num_points);
    if(m_ddc.size()) m_ddc.resize(num_components, num_points);
}

auto ChemicalField::val() -> Vector&
{
    return m_val;
}

auto ChemicalField::val() const -> const Vector&
{
    return m_val;
}

auto ChemicalField::val(Index i, double val) -> void
{
    m_val[i] = val;
}

auto ChemicalField::ddt(Index i, double ddt) -> void
{
    if(m_ddt.size()) m_ddt[i] = ddt;
}

auto ChemicalField::ddp(Index i, double ddp) -> void
{
    if(m_ddp.size()) m_ddp[i] = ddp;
}

auto ChemicalField::ddc(Index i, const Vector& ddc) -> void
{
    if(m_ddc.size()) m_ddc.col(i) = ddc;
}

auto ChemicalField::ddt(bool active) -> void
{
    if(active) m_ddt.resize(1);
    else m_ddt.resize(0);
}

auto ChemicalField::ddt() -> Vector&
{
    return m_ddt;
}

auto ChemicalField::ddt() const -> const Vector&
{
    return m_ddt;
}

auto ChemicalField::ddp(bool active) -> void
{
    if(active) m_ddp.resize(1);
    else m_ddp.resize(0);
}

auto ChemicalField::ddp() -> Vector&
{
    return m_ddp;
}

auto ChemicalField::ddp() const -> const Vector&
{
    return m_ddp;
}

auto ChemicalField::ddc(bool active) -> void
{
    if(active) m_ddc.resize(1, 1);
    else m_ddc.resize(0, 0);
}

auto ChemicalField::ddc() -> Matrix&
{
    return m_ddc;
}

auto ChemicalField::ddc() const -> const Matrix&
{
    return m_ddc;
}

}  // namespace Reaktoro
