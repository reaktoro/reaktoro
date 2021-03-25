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
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionThermoModelPhreeqcLgK class", "[ReactionThermoModelPhreeqcLgK]")
{
    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto A1 = 1.0;
    const auto A2 = 2.0;
    const auto A3 = 3.0;
    const auto A4 = 4.0;
    const auto A5 = 5.0;
    const auto A6 = 6.0;
    const auto Pr = 7.0;

    const auto model = ReactionThermoModelPhreeqcLgK({A1, A2, A3, A4, A5, A6, Pr});

    const auto R = universalGasConstant;

    const auto lgK   = A1 + A2*T + A3/T + A4*log10(T) + A5/(T*T) + A6*T*T;
    const auto lgK_T = A2 - A3/(T*T) + A4/(ln10*T) - 2*A5/(T*T*T) + 2*A6*T;

    const auto lnK   = ln10 * lgK;
    const auto lnK_T = ln10 * lgK_T;

    const auto dE = dV0 * (P - Pr);

    const auto dG0x = -R*T*lnK + dE;     // expected dG0 at (T, P)
    const auto dH0x =  R*T*T*lnK_T + dE; // expected dH0 at (T, P)

    ReactionThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0 == Approx(dG0x) );
    CHECK( rprops.dH0 == Approx(dH0x) );
}
