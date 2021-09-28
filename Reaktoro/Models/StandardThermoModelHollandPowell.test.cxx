// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelHollandPowell class", "[StandardThermoModelHollandPowell]")
{
    // Check Oelkers et al. (1995), page 1551, table for Quartz computed standard Gibbs energy.
    SECTION("testing standard thermodynamic properties for Quartz")
    {
        // Parameters for Quartz from SUPCRTBL (converted to SI units)
        StandardThermoModelParamsHollandPowell params;
        params.Gf       = -856280.0;
        params.Hf       = -910700.0;
        params.Sr       =  41.43;
        params.Vr       =  2.269e-05;
        params.a        =  92.9;
        params.b        = -0.000642;
        params.c        = -714900.0;
        params.d        = -716.1;
        params.alpha0   =  0.0;
        params.kappa0   =  73000000000.0;
        params.kappa0p  =  6.0;
        params.kappa0pp = -8.2e-11;
        params.numatoms =  3.0;

        auto model = StandardThermoModelHollandPowell(params);

        //======================================================================
        // Test method Model::operator()(T, P)
        //======================================================================

        WHEN("temperature is 75 C and pressure is 500 bar")
        {
            const auto T = 75.0 + 273.15; // in K
            const auto P = 500.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-204.923) ); // compare with -204.91 kcal/mol from Oelkers et al. (1995) p. 1551
            CHECK( props.H0  == Approx(-907269) );
            CHECK( props.V0  == Approx(2.26745e-05) );
            CHECK( props.Cp0 == Approx(48.3997) );
        }

        WHEN("temperature is 1000 C and pressure is 5000 bar")
        {
            const auto T = 1000.0 + 273.15; // in K
            const auto P = 5000.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-223.483) ); // compare with -223.91 kcal/mol from Oelkers et al. (1995) p. 1551

            CHECK( props.H0  == Approx(-837517) );
            CHECK( props.V0  == Approx(2.25382e-05) );
            CHECK( props.Cp0 == Approx(71.5722) );
        }

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().size() == 12 );
        CHECK( model.params()[0]  == params.Gf );
        CHECK( model.params()[1]  == params.Hf );
        CHECK( model.params()[2]  == params.Sr );
        CHECK( model.params()[3]  == params.Vr );
        CHECK( model.params()[4]  == params.a );
        CHECK( model.params()[5]  == params.b );
        CHECK( model.params()[6]  == params.c );
        CHECK( model.params()[7]  == params.d );
        CHECK( model.params()[8]  == params.alpha0 );
        CHECK( model.params()[9]  == params.kappa0 );
        CHECK( model.params()[10] == params.kappa0p );
        CHECK( model.params()[11] == params.kappa0pp );

        //======================================================================
        // Test method Model::serialize()
        //======================================================================

        yaml node;

        node = model.serialize();
        CHECK( double(node.at("HollandPowell").at("Gf"))       == params.Gf );
        CHECK( double(node.at("HollandPowell").at("Hf"))       == params.Hf );
        CHECK( double(node.at("HollandPowell").at("Sr"))       == params.Sr );
        CHECK( double(node.at("HollandPowell").at("Vr"))       == params.Vr );
        CHECK( double(node.at("HollandPowell").at("a"))        == params.a );
        CHECK( double(node.at("HollandPowell").at("b"))        == params.b );
        CHECK( double(node.at("HollandPowell").at("c"))        == params.c );
        CHECK( double(node.at("HollandPowell").at("d"))        == params.d );
        CHECK( double(node.at("HollandPowell").at("alpha0"))   == params.alpha0 );
        CHECK( double(node.at("HollandPowell").at("kappa0"))   == params.kappa0 );
        CHECK( double(node.at("HollandPowell").at("kappa0p"))  == params.kappa0p );
        CHECK( double(node.at("HollandPowell").at("kappa0pp")) == params.kappa0pp );

        params.Gf = 1234.0; // change value of Param object and check if new serialize call reflects this change

        node = model.serialize();
        CHECK( double(node.at("HollandPowell").at("Gf")) == 1234.0 );
    }

    // Check Oelkers et al. (1995), page 1544, table for Albite computed standard Gibbs energy.
    SECTION("testing standard thermodynamic properties for Albite")
    {
        // Parameters for Albite from SUPCRTBL (converted to SI units)
        StandardThermoModelParamsHollandPowell params;
        params.Gf       = -3712100.0;
        params.Hf       = -3935490.0;
        params.Sr       =  207.4;
        params.Vr       =  0.00010067;
        params.a        =  452.0;
        params.b        = -0.013364;
        params.c        = -1275900.0;
        params.d        = -3953.6;
        params.alpha0   =  2.36e-05;
        params.kappa0   =  54100000000.0;
        params.kappa0p  =  5.91;
        params.kappa0pp = -1.09e-10;
        params.numatoms =  13.0;

        auto model = StandardThermoModelHollandPowell(params);

        WHEN("temperature is 75 C and pressure is 500 bar")
        {
            const auto T = 75.0 + 273.15; // in K
            const auto P = 500.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-888.689) ); // compare with -890.42 kcal/mol from Oelkers et al. (1995) p. 1551
            CHECK( props.H0  == Approx(-3.91969e+06) );
            CHECK( props.V0  == Approx(0.000100699) );
            CHECK( props.Cp0 == Approx(224.931) );
        }

        WHEN("temperature is 1000 C and pressure is 5000 bar")
        {
            const auto T = 1000.0 + 273.15; // in K
            const auto P = 5000.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-977.639) ); // compare with -981.44 kcal/mol from Oelkers et al. (1995) p. 1551
            CHECK( props.H0  == Approx(-3.60235e+06) );
            CHECK( props.V0  == Approx(0.000102574) );
            CHECK( props.Cp0 == Approx(323.395) );
        }
    }

    // Check Oelkers et al. (1995), page 1553, table for CO2(g) computed standard Gibbs energy.
    SECTION("testing standard thermodynamic properties for CO2(g)")
    {
        // Parameters for CO2(g) from SUPCRTBL (converted to SI units)
        StandardThermoModelParamsHollandPowell params;
        params.Gf = -394350.0;
        params.Hf = -393510.0;
        params.Sr =  213.7;
        params.Vr =  0.0;
        params.a  =  87.8;
        params.b  = -0.002644;
        params.c  =  706400.0;
        params.d  = -998.9;

        auto model = StandardThermoModelHollandPowell(params);

        WHEN("temperature is 75 C and pressure is 500 bar")
        {
            const auto T = 75.0 + 273.15; // in K
            const auto P = 500.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-96.8416) ); // compare with -96.84 kcal/mol from Oelkers et al. (1995) p. 1553
            CHECK( props.H0  == Approx(-391603) );
            CHECK( props.V0  == Approx(0.0)     );
            CHECK( props.Cp0 == Approx(39.1723) );
        }

        WHEN("temperature is 1000 C and pressure is 5000 bar")
        {
            const auto T = 1000.0 + 273.15; // in K
            const auto P = 5000.0 * 1e5;  // in Pa

            StandardThermoProps props;
            props = model(T, P);

            CHECK( props.G0/4184  == Approx(-153.406) ); // compare with -153.42 kcal/mol from Oelkers et al. (1995) p. 1553
            CHECK( props.H0  == Approx(-344904) );
            CHECK( props.V0  == Approx(0.0)     );
            CHECK( props.Cp0 == Approx(56.8745) );
        }
    }
}
