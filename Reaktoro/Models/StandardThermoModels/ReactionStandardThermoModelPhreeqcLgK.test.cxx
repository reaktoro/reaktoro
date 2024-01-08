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
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionStandardThermoModelPhreeqcLgK class", "[ReactionStandardThermoModelPhreeqcLgK]")
{
    Param A1 = 1.0;
    Param A2 = 2.0;
    Param A3 = 3.0;
    Param A4 = 4.0;
    Param A5 = 5.0;
    Param A6 = 6.0;
    real  Pr = 7.0;

    const auto model = ReactionStandardThermoModelPhreeqcLgK({A1, A2, A3, A4, A5, A6, Pr});

    //======================================================================
    // Test method Model::operator()(T, P, dV0)
    //======================================================================

    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto R = universalGasConstant;

    const auto lgK    = A1 + A2*T + A3/T + A4*log10(T) + A5/(T*T) + A6*T*T;
    const auto lgK_T  = A2 - A3/(T*T) + A4/(ln10*T) - 2*A5/(T*T*T) + 2*A6*T;
    const auto lgK_TT = 2*A3/(T*T*T) - A4/(ln10*T*T) + 6*A5/(T*T*T*T) + 2*A6;

    const auto lnK    = ln10 * lgK;
    const auto lnK_T  = ln10 * lgK_T;
    const auto lnK_TT = ln10 * lgK_TT;

    const auto dE = dV0 * (P - Pr);

    const auto dG0x  = -R*T*lnK + dE;             // expected dG0 at (T, P)
    const auto dH0x  =  R*T*T*lnK_T + dE;         // expected dH0 at (T, P)
    const auto dCp0x =  2*R*T*lnK_T + R*T*T*lnK_TT; // expected dCp0 at (T, P)

    ReactionStandardThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0  == Approx(dG0x)  );
    CHECK( rprops.dH0  == Approx(dH0x)  );
    CHECK( rprops.dCp0 == Approx(dCp0x) );

    //======================================================================
    // Test method Model::params()
    //======================================================================

    CHECK( model.params().isDict() );
    CHECK( model.params().at("PhreeqcLgK").at("A1").asFloat() == A1 );
    CHECK( model.params().at("PhreeqcLgK").at("A2").asFloat() == A2 );
    CHECK( model.params().at("PhreeqcLgK").at("A3").asFloat() == A3 );
    CHECK( model.params().at("PhreeqcLgK").at("A4").asFloat() == A4 );
    CHECK( model.params().at("PhreeqcLgK").at("A5").asFloat() == A5 );
    CHECK( model.params().at("PhreeqcLgK").at("A6").asFloat() == A6 );
    CHECK( model.params().at("PhreeqcLgK").at("Pr").asFloat() == Pr );
}
