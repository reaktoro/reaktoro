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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelExtendedUNIQUAC.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelExtendedUNIQUAC class", "[StandardThermoModelExtendedUNIQUAC]")
{
    StandardThermoModelParamsExtendedUNIQUAC params;
    const auto Gr    = params.Gr    = 1.0;
    const auto Hr    = params.Hr    = 2.0;
    const auto Sr    = params.Sr    = 3.0;
    const auto Vr    = params.Vr    = 4.0;
    const auto Cp    = params.Cp    = 5.0;
    const auto a     = params.a     = 6.0;
    const auto b     = params.b     = 7.0;
    const auto c     = params.c     = 8.0;
    const auto alpha = params.alpha = 9.0;
    const auto beta  = params.beta  = 10.0;
    const auto Theta = params.Theta = 11.0;

    auto model = StandardThermoModelExtendedUNIQUAC(params);

    const auto Tr = 298.15;
    const auto Pr = 1.0;
    const auto KJ_TO_J = 1.0e3;
    const auto R = universalGasConstant;

    auto T = 360.0;
    auto P = 10.0e5;
    auto Pb = P * 1e-5;

    auto props = model(T, P);

    CHECK(props.Cp0 == Approx(a + b*T + c/(T - Theta)));
    CHECK(props.H0  == Approx(Hr*KJ_TO_J + a*(T - Tr) + 0.5*b*(T*T - Tr*Tr) + c*log((T - Theta)/(Tr - Theta))));
    CHECK(props.G0  == Approx(Gr*KJ_TO_J*T/Tr + Hr*KJ_TO_J*(1 - T/Tr) - a*T*(log(T/Tr) + Tr/T - 1) - 0.5*b*(T - Tr)*(T - Tr) - c*T/Theta*((T - Theta)/T*log((T - Theta)/(Tr - Theta)) - log(T/Tr)) + R*T*(alpha + beta*(Pb - Pr))*(Pb - Pr)));

    //======================================================================
    // Test method Model::params()
    //======================================================================

    CHECK( model.params().isDict() );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Gr").asFloat()    == params.Gr );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Hr").asFloat()    == params.Hr );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Sr").asFloat()    == params.Sr );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Vr").asFloat()    == params.Vr );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Cp").asFloat()    == params.Cp );
    CHECK( model.params().at("ExtendedUNIQUAC").at("a").asFloat()     == params.a );
    CHECK( model.params().at("ExtendedUNIQUAC").at("b").asFloat()     == params.b );
    CHECK( model.params().at("ExtendedUNIQUAC").at("c").asFloat()     == params.c );
    CHECK( model.params().at("ExtendedUNIQUAC").at("alpha").asFloat() == params.alpha );
    CHECK( model.params().at("ExtendedUNIQUAC").at("beta").asFloat()  == params.beta );
    CHECK( model.params().at("ExtendedUNIQUAC").at("Theta").asFloat() == params.Theta );
}
