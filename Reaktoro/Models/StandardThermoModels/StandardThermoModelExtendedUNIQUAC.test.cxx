// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

    CHECK( model.params().size() == 11 );
    CHECK( model.params()[0]  == params.Gr );
    CHECK( model.params()[1]  == params.Hr );
    CHECK( model.params()[2]  == params.Sr );
    CHECK( model.params()[3]  == params.Vr );
    CHECK( model.params()[4]  == params.Cp );
    CHECK( model.params()[5]  == params.a );
    CHECK( model.params()[6]  == params.b );
    CHECK( model.params()[7]  == params.c );
    CHECK( model.params()[8]  == params.alpha );
    CHECK( model.params()[9]  == params.beta );
    CHECK( model.params()[10] == params.Theta );

    //======================================================================
    // Test method Model::serialize()
    //======================================================================

    Data node;

    node = model.serialize();
    CHECK( double(node.at("ExtendedUNIQUAC").at("Gr"))    == params.Gr );
    CHECK( double(node.at("ExtendedUNIQUAC").at("Hr"))    == params.Hr );
    CHECK( double(node.at("ExtendedUNIQUAC").at("Sr"))    == params.Sr );
    CHECK( double(node.at("ExtendedUNIQUAC").at("Vr"))    == params.Vr );
    CHECK( double(node.at("ExtendedUNIQUAC").at("Cp"))    == params.Cp );
    CHECK( double(node.at("ExtendedUNIQUAC").at("a"))     == params.a );
    CHECK( double(node.at("ExtendedUNIQUAC").at("b"))     == params.b );
    CHECK( double(node.at("ExtendedUNIQUAC").at("c"))     == params.c );
    CHECK( double(node.at("ExtendedUNIQUAC").at("alpha")) == params.alpha );
    CHECK( double(node.at("ExtendedUNIQUAC").at("beta"))  == params.beta );
    CHECK( double(node.at("ExtendedUNIQUAC").at("Theta")) == params.Theta );

    params.Gr = 1234.0; // change value of Param object and check if new serialize call reflects this change

    node = model.serialize();

    CHECK( double(node.at("ExtendedUNIQUAC").at("Gr")) == 1234.0 );
}
