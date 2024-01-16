// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "WaterThermoProps.hpp"

#define WATER_THERMO_PROPS_ADD_SUB_IMPL(op) \
    T   op x.T;                             \
    V   op x.V;                             \
    S   op x.S;                             \
    A   op x.A;                             \
    U   op x.U;                             \
    H   op x.H;                             \
    G   op x.G;                             \
    Cv  op x.Cv;                            \
    Cp  op x.Cp;                            \
    D   op x.D;                             \
    DT  op x.DT;                            \
    DP  op x.DP;                            \
    DTT op x.DTT;                           \
    DTP op x.DTP;                           \
    DPP op x.DPP;                           \
    P   op x.P;                             \
    PT  op x.PT;                            \
    PD  op x.PD;                            \
    PTT op x.PTT;                           \
    PTD op x.PTD;                           \
    PDD op x.PDD;                           \
    return *this;

#define WATER_THERMO_PROPS_MUL_DIV_IMPL(op) \
    T   op x;                               \
    V   op x;                               \
    S   op x;                               \
    A   op x;                               \
    U   op x;                               \
    H   op x;                               \
    G   op x;                               \
    Cv  op x;                               \
    Cp  op x;                               \
    D   op x;                               \
    DT  op x;                               \
    DP  op x;                               \
    DTT op x;                               \
    DTP op x;                               \
    DPP op x;                               \
    P   op x;                               \
    PT  op x;                               \
    PD  op x;                               \
    PTT op x;                               \
    PTD op x;                               \
    PDD op x;                               \
    return *this;

namespace Reaktoro {

auto WaterThermoProps::operator+=(WaterThermoProps const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_ADD_SUB_IMPL(+=)
}

auto WaterThermoProps::operator-=(WaterThermoProps const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_ADD_SUB_IMPL(-=)
}

auto WaterThermoProps::operator*=(double const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_MUL_DIV_IMPL(*=)
}

auto WaterThermoProps::operator*=(real const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_MUL_DIV_IMPL(*=)
}

auto WaterThermoProps::operator/=(double const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_MUL_DIV_IMPL(/=)
}

auto WaterThermoProps::operator/=(real const& x) -> WaterThermoProps&
{
    WATER_THERMO_PROPS_MUL_DIV_IMPL(/=)
}

auto operator+(WaterThermoProps const& r) -> WaterThermoProps
{
    return r;
}

auto operator+(WaterThermoProps&& r) -> WaterThermoProps&&
{
    return std::move(r);
}

auto operator-(WaterThermoProps const& r) -> WaterThermoProps
{
    auto result(r);
    return result *= -1.0;
}

auto operator-(WaterThermoProps&& r) -> WaterThermoProps&&
{
    r *= -1.0;
    return std::move(r);
}

auto operator+(WaterThermoProps const& l, WaterThermoProps const& r) -> WaterThermoProps
{
    auto result(l);
    return result += r;
}

auto operator+(WaterThermoProps&& l, WaterThermoProps const& r) -> WaterThermoProps&&
{
    l += r;
    return std::move(l);
}

auto operator+(WaterThermoProps const& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    r += l;
    return std::move(r);
}

auto operator+(WaterThermoProps&& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    l += r;
    return std::move(l);
}

auto operator-(WaterThermoProps const& l, WaterThermoProps const& r) -> WaterThermoProps
{
    auto result(l);
    return result -= r;
}

auto operator-(WaterThermoProps&& l, WaterThermoProps const& r) -> WaterThermoProps&&
{
    l -= r;
    return std::move(l);
}

auto operator-(WaterThermoProps const& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    r *= -1.0;
    r += l;
    return std::move(r);
}

auto operator-(WaterThermoProps&& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    l -= r;
    return std::move(l);
}

auto operator*(double const& l, WaterThermoProps const& r) -> WaterThermoProps
{
    auto result(r);
    return result *= l;
}

auto operator*(double const& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    r *= l;
    return std::move(r);
}

auto operator*(WaterThermoProps const& l, double const& r) -> WaterThermoProps
{
    return r * l;
}

auto operator*(WaterThermoProps&& l, double const& r) -> WaterThermoProps&&
{
    l *= r;
    return std::move(l);
}

auto operator*(real const& l, WaterThermoProps const& r) -> WaterThermoProps
{
    auto result(r);
    return result *= l;
}

auto operator*(real const& l, WaterThermoProps&& r) -> WaterThermoProps&&
{
    r *= l;
    return std::move(r);
}

auto operator*(WaterThermoProps const& l, real const& r) -> WaterThermoProps
{
    return r * l;
}

auto operator*(WaterThermoProps&& l, real const& r) -> WaterThermoProps&&
{
    l *= r;
    return std::move(l);
}

auto operator/(WaterThermoProps const& l, double const& r) -> WaterThermoProps
{
    auto result(l);
    return result /= r;
}

auto operator/(WaterThermoProps&& l, double const& r) -> WaterThermoProps&&
{
    l /= r;
    return std::move(l);
}

auto operator/(WaterThermoProps const& l, real const& r) -> WaterThermoProps
{
    auto result(l);
    return result /= r;
}

auto operator/(WaterThermoProps&& l, real const& r) -> WaterThermoProps&&
{
    l /= r;
    return std::move(l);
}

} // namespace Reaktoro
