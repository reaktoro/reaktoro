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
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionStandardThermoModelVantHoff class", "[ReactionStandardThermoModelVantHoff]")
{
    real lgKr = 1.0;
    real dHr  = 2.0;
    real Tr   = 3.0;
    real Pr   = 4.0;

    const auto model = ReactionStandardThermoModelVantHoff({lgKr, dHr, Tr, Pr});

    //======================================================================
    // Test method Model::operator()(T, P, dV0)
    //======================================================================

    const auto T = 5.0;
    const auto P = 7.0;
    const auto dV0 = 9.0;

    const auto R = universalGasConstant;

    const auto lnKr = lgKr * ln10;

    const auto lnK = lnKr - dHr/R * (1.0/T - 1.0/Tr);

    const auto dE = dV0*(P - Pr);

    const auto dG0x = -R*T*lnK + dE; // expected dG0 at (T, P)
    const auto dH0x = dHr + dE; // expected dH0 at (T, P)

    ReactionStandardThermoProps rprops = model({T, P, dV0});

    CHECK( rprops.dG0  == Approx(dG0x) );
    CHECK( rprops.dH0  == Approx(dH0x) );
    CHECK( rprops.dCp0 == 0.0 );

    //======================================================================
    // Test method Model::params()
    //======================================================================

    CHECK( model.params().isDict() );
    CHECK( model.params().at("VantHoff").at("lgKr").asFloat() == lgKr );
    CHECK( model.params().at("VantHoff").at("dHr").asFloat()  == dHr );
    CHECK( model.params().at("VantHoff").at("Tr").asFloat()   == Tr );
    CHECK( model.params().at("VantHoff").at("Pr").asFloat()   == Pr );
}
