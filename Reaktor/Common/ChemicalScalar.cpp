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
#include <Reaktor/Common/ThermoScalar.hpp>

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

ChemicalScalar::ChemicalScalar(const ChemicalVectorRowConst& row)
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

auto ChemicalScalar::operator=(const ChemicalVectorRowConst& row) -> ChemicalScalar&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    ddn = row.ddn;
    return *this;
}

auto ChemicalScalar::operator+=(const ChemicalScalar& other) -> ChemicalScalar&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    ddn += other.ddn;
    return *this;
}

auto ChemicalScalar::operator+=(const ThermoScalar& other) -> ChemicalScalar&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    return *this;
}

auto ChemicalScalar::operator-=(const ChemicalScalar& other) -> ChemicalScalar&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    ddn -= other.ddn;
    return *this;
}

auto ChemicalScalar::operator-=(const ThermoScalar& other) -> ChemicalScalar&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    return *this;
}

auto ChemicalScalar::operator*=(double scalar) -> ChemicalScalar&
{
    val *= scalar;
    ddt *= scalar;
    ddp *= scalar;
    ddn *= scalar;
    return *this;
}

auto ChemicalScalar::operator/=(double scalar) -> ChemicalScalar&
{
    *this *= 1.0/scalar;
    return *this;
}

auto operator==(const ChemicalScalar& l, const ChemicalScalar& r) -> bool
{
    return l.val == r.val and
           l.ddt == r.ddt and
           l.ddp == r.ddp and
           l.ddn == r.ddn;
}

auto operator+(const ChemicalScalar& l) -> ChemicalScalar
{
    return l;
}

auto operator-(const ChemicalScalar& l) -> ChemicalScalar
{
    return -1.0 * l;
}

auto operator+(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = l;
    res += r;
    return res;
}

auto operator+(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = l;
    res += r;
    return res;
}

auto operator+(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = r;
    res += l;
    return res;
}

auto operator-(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = l;
    res -= r;
    return res;
}

auto operator-(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = l;
    res -= r;
    return res;
}

auto operator-(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = r;
    res -= l;
    return res;
}

auto operator*(double scalar, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res = r;
    res *= scalar;
    return res;
}

auto operator*(const ChemicalScalar& l, double scalar) -> ChemicalScalar
{
    return scalar * l;
}

auto operator*(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    ChemicalScalar res;
    res.val = l.val * r.val;
    res.ddt = l.ddt * r.val + l.val * r.ddt;
    res.ddp = l.ddp * r.val + l.val * r.ddp;
    res.ddn = l.ddn * r.val + l.val * r.ddn;
    return res;
}

auto operator/(double scalar, const ChemicalScalar& r) -> ChemicalScalar
{
    const double factor = -scalar/(r.val * r.val);
    ChemicalScalar res;
    res.val = scalar/r.val;
    res.ddt = factor * r.ddt;
    res.ddp = factor * r.ddp;
    res.ddn = factor * r.ddn;
    return res;
}

auto operator/(const ChemicalScalar& l, double scalar) -> ChemicalScalar
{
    return (1.0/scalar) * l;
}

auto operator/(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar
{
    const double factor = 1.0/(r.val * r.val);
    ChemicalScalar res;
    res.val = l.val / r.val;
    res.ddt = (r.val * l.ddt - l.val * r.ddt) * factor;
    res.ddp = (r.val * l.ddp - l.val * r.ddp) * factor;
    res.ddn = (r.val * l.ddn - l.val * r.ddn) * factor;
    return res;
}

auto pow(const ChemicalScalar& a, double power) -> ChemicalScalar
{
    ChemicalScalar b;
    b.val = std::pow(a.val, power);
    b.ddt = power * b.val * a.ddt/a.val;
    b.ddp = power * b.val * a.ddp/a.val;
    b.ddn = power * b.val * a.ddn/a.val;
    return b;
}

auto exp(const ChemicalScalar& a) -> ChemicalScalar
{
    ChemicalScalar b;
    b.val = std::exp(a.val);
    b.ddt = b.val * a.ddt;
    b.ddp = b.val * a.ddp;
    b.ddn = b.val * a.ddn;
    return b;
}

auto log(const ChemicalScalar& a) -> ChemicalScalar
{
    const double factor = 1.0/a.val;
    ChemicalScalar b;
    b.val = std::log(a.val);
    b.ddt = factor * a.ddt;
    b.ddp = factor * a.ddp;
    b.ddn = factor * a.ddn;
    return b;
}

} // namespace Reaktor
