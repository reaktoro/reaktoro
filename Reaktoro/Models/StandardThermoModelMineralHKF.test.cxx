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
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelMineralHKF class", "[StandardThermoModelMineralHKF]")
{
    double T = 0.0;
    double P = 0.0;

    // Check Oelkers et al. (1995), page 1544, table for Albite-low.
    SECTION("testing standard thermodynamic properties for Albite-low")
    {
        // Parameters for Albite,low from slop98.dat (converted to SI units)
        StandardThermoModelParamsMineralHKF params;
        params.Gf   = -3708312.7;
        params.Hf   = -3931621.1;
        params.Sr   =  207.14984;
        params.Vr   =  0.00010007;
        params.ntr  =  0;
        params.a    =  {258.1528};
        params.b    =  {0.0581576};
        params.c    =  {-6280184.0};
        params.Tmax =  2500.0;

        auto model = StandardThermoModelMineralHKF(params);

        //======================================================================
        // Test method Model::operator()(T, P)
        //======================================================================

        StandardThermoProps props;

        T = 75.0 + 273.15; // 75 degC (in K)
        P = 500.0 * 1e5;  // 500 bar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-3.71452e+06) ); // converted to J/mol from -890.42 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.91581e+06) );
        CHECK( props.V0  == Approx(0.00010007)   );
        CHECK( props.Cp0 == Approx(226.587)      );
        CHECK( props.Cv0 == Approx(226.587)      );

        T = 1000.0 + 273.15; // 1000 degC (in K)
        P = 5000.0 * 1e5;  // 5kbar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-4.08694e+06) ); // converted to J/mol from -979.43 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.60148e+06) );
        CHECK( props.V0  == Approx(0.00010007)   );
        CHECK( props.Cp0 == Approx(328.322)      );
        CHECK( props.Cv0 == Approx(328.322)      );

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().size() == 7 );
        CHECK( model.params()[0]  == params.Gf );
        CHECK( model.params()[1]  == params.Hf );
        CHECK( model.params()[2]  == params.Sr );
        CHECK( model.params()[3]  == params.Vr );
        CHECK( model.params()[4]  == params.a[0] );
        CHECK( model.params()[5]  == params.b[0] );
        CHECK( model.params()[6]  == params.c[0] );

        //======================================================================
        // Test method Model::serialize()
        //======================================================================

        yaml node;

        node = model.serialize();
        CHECK( double(node.at("MineralHKF").at("Gf"))   == params.Gf );
        CHECK( double(node.at("MineralHKF").at("Hf"))   == params.Hf );
        CHECK( double(node.at("MineralHKF").at("Sr"))   == params.Sr );
        CHECK( double(node.at("MineralHKF").at("Vr"))   == params.Vr );
        CHECK( double(node.at("MineralHKF").at("a")[0]) == params.a[0] );
        CHECK( double(node.at("MineralHKF").at("b")[0]) == params.b[0] );
        CHECK( double(node.at("MineralHKF").at("c")[0]) == params.c[0] );

        params.Gf = 1234.0; // change value of Param object and check if new serialize call reflects this change

        node = model.serialize();
        CHECK( double(node.at("MineralHKF").at("Gf")) == 1234.0 );
    }

    // Check Oelkers et al. (1995), page 1544, table for Albite.
    SECTION("testing standard thermodynamic properties for Albite")
    {
        // Parameters for Albite,low from slop98.dat (converted to SI units)
        StandardThermoModelParamsMineralHKF params;
        params.Gf = -3708312.7;
        params.Hf = -3931621.1;
        params.Sr = 207.14984;
        params.Vr = 0.00010025;
        params.ntr = 1;
        params.a = {258.1528, 342.58592};
        params.b = {0.0581576, 0.014869936};
        params.c = {-6280184.0, -20984434.0};
        params.Ttr = {473.0};
        params.Htr = {0.0};
        params.Vtr = {0.0};
        params.dPdTtr = {0.0};

        auto model = StandardThermoModelMineralHKF(params);

        //======================================================================
        // Test method Model::operator()(T, P)
        //======================================================================

        StandardThermoProps props;

        T = 75.0 + 273.15; // 75 degC (in K)
        P = 500.0 * 1e5;  // 500 bar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-3.71451e+06) ); // converted to J/mol from -890.42 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.9158e+06)  );
        CHECK( props.V0  == Approx(0.00010025)   );
        CHECK( props.Cp0 == Approx(226.587)      );
        CHECK( props.Cv0 == Approx(226.587)      );

        T = 1000.0 + 273.15; // 1000 degC (in K)
        P = 5000.0 * 1e5;  // 5kbar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-4.09534e+06) ); // converted to J/mol from -979.43 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.58361e+06) );
        CHECK( props.V0  == Approx(0.00010025)   );
        CHECK( props.Cp0 == Approx(348.572)      );
        CHECK( props.Cv0 == Approx(348.572)      );

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().size() == 14 );
        CHECK( model.params()[0]  == params.Gf );
        CHECK( model.params()[1]  == params.Hf );
        CHECK( model.params()[2]  == params.Sr );
        CHECK( model.params()[3]  == params.Vr );
        CHECK( model.params()[4]  == params.a[0] );
        CHECK( model.params()[5]  == params.a[1] );
        CHECK( model.params()[6]  == params.b[0] );
        CHECK( model.params()[7]  == params.b[1] );
        CHECK( model.params()[8]  == params.c[0] );
        CHECK( model.params()[9]  == params.c[1] );
        CHECK( model.params()[10] == params.Ttr[0] );
        CHECK( model.params()[11] == params.Htr[0] );
        CHECK( model.params()[12] == params.Vtr[0] );
        CHECK( model.params()[13] == params.dPdTtr[0] );

        //======================================================================
        // Test method Model::serialize()
        //======================================================================

        yaml node;

        node = model.serialize();
        CHECK( double(node.at("MineralHKF").at("Gf"))        == params.Gf );
        CHECK( double(node.at("MineralHKF").at("Hf"))        == params.Hf );
        CHECK( double(node.at("MineralHKF").at("Sr"))        == params.Sr );
        CHECK( double(node.at("MineralHKF").at("Vr"))        == params.Vr );
        CHECK( double(node.at("MineralHKF").at("a")[0])      == params.a[0] );
        CHECK( double(node.at("MineralHKF").at("a")[1])      == params.a[1] );
        CHECK( double(node.at("MineralHKF").at("b")[0])      == params.b[0] );
        CHECK( double(node.at("MineralHKF").at("b")[1])      == params.b[1] );
        CHECK( double(node.at("MineralHKF").at("c")[0])      == params.c[0] );
        CHECK( double(node.at("MineralHKF").at("c")[1])      == params.c[1] );
        CHECK( double(node.at("MineralHKF").at("Ttr")[0])    == params.Ttr[0] );
        CHECK( double(node.at("MineralHKF").at("Htr")[0])    == params.Htr[0] );
        CHECK( double(node.at("MineralHKF").at("Vtr")[0])    == params.Vtr[0] );
        CHECK( double(node.at("MineralHKF").at("dPdTtr")[0]) == params.dPdTtr[0] );

        params.Gf = 1234.0; // change value of Param object and check if new serialize call reflects this change

        node = model.serialize();
        CHECK( double(node.at("MineralHKF").at("Gf")) == 1234.0 );
    }
}
