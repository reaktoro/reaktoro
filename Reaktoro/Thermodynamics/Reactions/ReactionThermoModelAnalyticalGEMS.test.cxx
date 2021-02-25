// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelAnalyticalGEMS.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionThermoModelAnalyticalGEMS class", "[ReactionThermoModelAnalyticalGEMS]")
{
    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto A0 = 1.0;
    const auto A1 = 2.0;
    const auto A2 = 3.0;
    const auto A3 = 4.0;
    const auto A4 = 5.0;
    const auto A5 = 6.0;
    const auto A6 = 7.0;

    const auto model = ReactionThermoModelAnalyticalGEMS(A0, A1, A2, A3, A4, A5, A6);

    const auto R = universalGasConstant;

    const auto lgK   = A0 + A1*T + A2/T + A3*log(T) + A4/(T*T) + A5*(T*T) + A6/sqrt(T);
    const auto lgK_T = A1 - A2/(T*T) + A3/T - 2*A4/(T*T*T) + 2*A5*T - 0.5*A6/(T*sqrt(T));

    const auto lnK   = ln10 * lgK;
    const auto lnK_T = ln10 * lgK_T;

    const auto dG0x = -R*T*lnK;     // expected dG0 at (T, P)
    const auto dH0x =  R*T*T*lnK_T; // expected dH0 at (T, P)

    ReactionThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0 == Approx(dG0x) );
    CHECK( rprops.dH0 == Approx(dH0x) );
}
