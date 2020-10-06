// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Thermodynamics/Reactions/ReactionThermoModelVantHoff.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionThermoModelVantHoff class", "[ReactionThermoModelVantHoff]")
{
    const auto lgK0 = 1.0;
    const auto dH0  = 2.0;
    const auto Tref = 3.0;

    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto model = ReactionThermoModelVantHoff(lgK0, dH0, Tref);

    const auto rpropsfn = model;

    const auto R = universalGasConstant;

    const auto lnK0 = lgK0 * ln10;

    const auto lnK = lnK0 - dH0/R * (1.0/T - 1.0/Tref);

    const auto dG0x = -R*T*lnK; // expected dG0 at (T, P)
    const auto dH0x = dH0; // expected dH0 at (T, P)

    ReactionThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0 == Approx(dG0x) );
    CHECK( rprops.dH0 == Approx(dH0x) );
}
