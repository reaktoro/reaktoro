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
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelConstLgK.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionThermoModelConstLgK class", "[ReactionThermoModelConstLgK]")
{
    const auto lgKr = 1.0;
    const auto Pr   = 2.0;

    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto model = ReactionThermoModelConstLgK({lgKr, Pr});

    const auto rpropsfn = model;

    const auto R = universalGasConstant;

    const auto lnKr = lgKr * ln10;

    const auto dE = dV0 * (P - Pr);

    const auto dG0x = -R*T*lnKr + dE; // expected dG0 at (T, P)
    const auto dH0x = dE;             // expected dH0 at (T, P)

    ReactionThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0 == Approx(dG0x) );
    CHECK( rprops.dH0 == Approx(dH0x) );
}
