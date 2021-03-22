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
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelHKF.hpp>
using namespace Reaktoro;

//======================================================================
// NOTE:
//======================================================================
// The tests below use data from the following reference:
// Oelkers, E. H., Helgeson, H. C., Shock, E. L., Sverjensky, D. A., Johnson,
// J. W., & Pokrovskii, V. A. (1995). Summary of the Apparent Standard Partial
// Molal Gibbs Free Energies of Formation of Aqueous Species, Minerals, and
// Gases at Pressures 1 to 5000 Bars and Temperatures 25 to 1000 °C. Journal of
// Physical and Chemical Reference Data, 24(4), 1401.
// https://doi.org/10.1063/1.555976
//======================================================================

TEST_CASE("Testing StandardThermoModelHKF class", "[StandardThermoModelHKF]")
{
    const auto T = 75.0 + 273.15; // 75 degC (in K)
    const auto P = 1000.0 * 1e5;  // 1kbar (in Pa)

    // Check Oelkers et al. (1995), page 1438, table for CO2(aq).
    SECTION("testing standard thermodynamic properties for CO2(aq)")
    {
        // Parameters for CO2(aq) from slop98.dat (converted to SI units)
        StandardThermoModelParamsHKF params;
        params.Gf     = -385974.0;
        params.Hf     = -413797.6;
        params.Sr     =  117.5704;
        params.a1     =  2.6135774e-05;
        params.a2     =  3125.9082;
        params.a3     =  0.00011772102;
        params.a4     = -129197.74;
        params.c1     =  167.49598;
        params.c2     =  368208.74;
        params.wref   = -8368.0;
        params.charge =  0.0;
        params.formula = "CO2";

        auto model = StandardThermoModelHKF(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-389346.0)   ); // converted to J/mol from -93.06 kcal/mol as in table
        CHECK( props.H0  == Approx(-400560.0)   );
        CHECK( props.V0  == Approx(3.28667e-05) );
        CHECK( props.Cp0 == Approx(205.903)     );
    }

    // Check Oelkers et al. (1995), page 1438, table for CO3-2.
    SECTION("testing standard thermodynamic properties for CO3-2")
    {
        // Parameters for CO3-2 from slop98.dat (converted to SI units)
        StandardThermoModelParamsHKF params;
        params.Gf     = -527983.14;
        params.Hf     = -675234.84;
        params.Sr     = -49.9988;
        params.a1     =  1.1934442e-05;
        params.a2     = -1667.073;
        params.a3     =  0.00026837013;
        params.a4     = -109382.31;
        params.c1     = -13.89339;
        params.c2     = -719300.73;
        params.wref   =  1418961.8;
        params.charge = -2.0;
        params.formula = "CO3-2";

        auto model = StandardThermoModelHKF(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-524990.0)   ); // converted to J/mol from -125.48 kcal/mol as in table
        CHECK( props.H0  == Approx(-685294.0)   );
        CHECK( props.V0  == Approx(-2.3192e-06) );
        CHECK( props.Cp0 == Approx(-189.733)    );
    }

    // Check Oelkers et al. (1995), page 1463, table for H+.
    SECTION("testing standard thermodynamic properties for H+")
    {
        // Parameters for H+ from slop98.dat (converted to SI units)
        StandardThermoModelParamsHKF params;
        params.Gf     = 0.0;
        params.Hf     = 0.0;
        params.Sr     = 0.0;
        params.a1     = 0.0;
        params.a2     = 0.0;
        params.a3     = 0.0;
        params.a4     = 0.0;
        params.c1     = 0.0;
        params.c2     = 0.0;
        params.wref   = 0.0;
        params.charge = 1.0;
        params.formula = "H+";

        auto model = StandardThermoModelHKF(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(0.0) ); // converted to J/mol from 0.0 kcal/mol as in table
        CHECK( props.H0  == Approx(0.0) );
        CHECK( props.V0  == Approx(0.0) );
        CHECK( props.Cp0 == Approx(0.0) );
    }

    // Check Oelkers et al. (1995), page 1489, table for Mg+2.
    SECTION("testing standard thermodynamic properties for Mg+2")
    {
        // Parameters for Mg+2 from slop98.dat (converted to SI units)
        StandardThermoModelParamsHKF params;
        params.Gf     = -453984.92;
        params.Hf     = -465959.53;
        params.Sr     = -138.072;
        params.a1     = -3.4379928e-06;
        params.a2     = -3597.8216;
        params.a3     =  0.0003510376;
        params.a4     = -99997.6;
        params.c1     =  87.0272;
        params.c2     = -246521.28;
        params.wref   =  643164.48;
        params.charge =  2.0;
        params.formula = "Mg+2";

        auto model = StandardThermoModelHKF(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-449027.0)    ); // converted to J/mol from -107.32 kcal/mol as in table
        CHECK( props.H0  == Approx(-466977.0)    );
        CHECK( props.V0  == Approx(-1.70501e-05) );
        CHECK( props.Cp0 == Approx(10.2122)      );
    }
}
