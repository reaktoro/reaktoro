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

#include "ThermoScalar.hpp"

// Reaktor includes
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

ThermoScalar::ThermoScalar()
{}

ThermoScalar::ThermoScalar(double val, double ddt, double ddp)
: val(val), ddt(ddt), ddp(ddp)
{}

auto ThermoScalar::operator=(const ThermoVectorRow& row) -> ThermoScalar&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    return *this;
}

auto ThermoScalar::operator=(const ThermoVectorConstRow& row) -> ThermoScalar&
{
    val = row.val;
    ddt = row.ddt;
    ddp = row.ddp;
    return *this;
}

auto ThermoScalar::operator+=(const ThermoScalar& other) -> ThermoScalar&
{
    val += other.val;
    ddt += other.ddt;
    ddp += other.ddp;
    return *this;
}

auto ThermoScalar::operator-=(const ThermoScalar& other) -> ThermoScalar&
{
    val -= other.val;
    ddt -= other.ddt;
    ddp -= other.ddp;
    return *this;
}

auto ThermoScalar::operator*=(double scalar) -> ThermoScalar&
{
    val *= scalar;
    ddt *= scalar;
    ddp *= scalar;
    return *this;
}

auto ThermoScalar::operator/=(double scalar) -> ThermoScalar&
{
    *this *= 1.0/scalar;
    return *this;
}

auto operator==(const ThermoScalar& l, const ThermoScalar& r) -> bool
{
    return l.val == r.val and
           l.ddt == r.ddt and
           l.ddp == r.ddp;
}

auto operator+(const ThermoScalar& l) -> ThermoScalar
{
    return l;
}

auto operator-(const ThermoScalar& l) -> ThermoScalar
{
    return -1.0 * l;
}

auto operator+(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar
{
    ThermoScalar res = l;
    res += r;
    return res;
}

auto operator-(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar
{
    ThermoScalar res = l;
    res -= r;
    return res;
}

auto operator*(double scalar, const ThermoScalar& r) -> ThermoScalar
{
    ThermoScalar res = r;
    res *= scalar;
    return res;
}

auto operator*(const ThermoScalar& l, double scalar) -> ThermoScalar
{
    return scalar * l;
}

auto operator/(const ThermoScalar& l, double scalar) -> ThermoScalar
{
    return (1.0/scalar) * l;
}

} // namespace Reaktor
