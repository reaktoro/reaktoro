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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
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
        CHECK( props.VT0 == Approx(0.00000000)   );
        CHECK( props.VP0 == Approx(0.00000000)   );
        CHECK( props.Cp0 == Approx(226.587)      );

        T = 1000.0 + 273.15; // 1000 degC (in K)
        P = 5000.0 * 1e5;  // 5kbar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-4.08694e+06) ); // converted to J/mol from -979.43 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.60148e+06) );
        CHECK( props.V0  == Approx(0.00010007)   );
        CHECK( props.VT0 == Approx(0.00000000)   );
        CHECK( props.VP0 == Approx(0.00000000)   );
        CHECK( props.Cp0 == Approx(328.322)      );

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().isDict() );
        CHECK( model.params().at("MineralHKF").at("Gf").asFloat()   == params.Gf );
        CHECK( model.params().at("MineralHKF").at("Hf").asFloat()   == params.Hf );
        CHECK( model.params().at("MineralHKF").at("Sr").asFloat()   == params.Sr );
        CHECK( model.params().at("MineralHKF").at("Vr").asFloat()   == params.Vr );
        CHECK( model.params().at("MineralHKF").at("a")[0].asFloat() == params.a[0] );
        CHECK( model.params().at("MineralHKF").at("b")[0].asFloat() == params.b[0] );
        CHECK( model.params().at("MineralHKF").at("c")[0].asFloat() == params.c[0] );
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
        CHECK( props.VT0 == Approx(0.00000000)   );
        CHECK( props.VP0 == Approx(0.00000000)   );
        CHECK( props.Cp0 == Approx(226.587)      );

        T = 1000.0 + 273.15; // 1000 degC (in K)
        P = 5000.0 * 1e5;  // 5kbar (in Pa)
        props = model(T, P);

        CHECK( props.G0  == Approx(-4.09534e+06) ); // converted to J/mol from -979.43 kcal/mol as in table (there are some 0.3%-0.5% difference)
        CHECK( props.H0  == Approx(-3.58361e+06) );
        CHECK( props.V0  == Approx(0.00010025)   );
        CHECK( props.VT0 == Approx(0.00000000)   );
        CHECK( props.VP0 == Approx(0.00000000)   );
        CHECK( props.Cp0 == Approx(348.572)      );

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().isDict() );
        CHECK( model.params().at("MineralHKF").at("Gf").asFloat()        == params.Gf );
        CHECK( model.params().at("MineralHKF").at("Hf").asFloat()        == params.Hf );
        CHECK( model.params().at("MineralHKF").at("Sr").asFloat()        == params.Sr );
        CHECK( model.params().at("MineralHKF").at("Vr").asFloat()        == params.Vr );
        CHECK( model.params().at("MineralHKF").at("a")[0].asFloat()      == params.a[0] );
        CHECK( model.params().at("MineralHKF").at("a")[1].asFloat()      == params.a[1] );
        CHECK( model.params().at("MineralHKF").at("b")[0].asFloat()      == params.b[0] );
        CHECK( model.params().at("MineralHKF").at("b")[1].asFloat()      == params.b[1] );
        CHECK( model.params().at("MineralHKF").at("c")[0].asFloat()      == params.c[0] );
        CHECK( model.params().at("MineralHKF").at("c")[1].asFloat()      == params.c[1] );
        CHECK( model.params().at("MineralHKF").at("Ttr")[0].asFloat()    == params.Ttr[0] );
        CHECK( model.params().at("MineralHKF").at("Htr")[0].asFloat()    == params.Htr[0] );
        CHECK( model.params().at("MineralHKF").at("Vtr")[0].asFloat()    == params.Vtr[0] );
        CHECK( model.params().at("MineralHKF").at("dPdTtr")[0].asFloat() == params.dPdTtr[0] );
    }
}
