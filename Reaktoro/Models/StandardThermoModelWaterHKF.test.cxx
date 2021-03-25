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
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelWaterHKF class", "[StandardThermoModelWaterHKF]")
{
    const auto T = 75.0 + 273.15; // 75 degC (in K)
    const auto P = 1000.0 * 1e5;  // 1kbar (in Pa)

    // Check Oelkers et al. (1995), page 1542, table for H2O.
    SECTION("testing standard thermodynamic properties for H2O")
    {
        // Parameters for H2O from slop98.dat (converted to SI units)
        StandardThermoModelParamsWaterHKF params;
        params.Ttr = 273.16;
        params.Str = 63.312288;
        params.Gtr = -235517.36;
        params.Htr = -287721.128;

        auto model = StandardThermoModelWaterHKF(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-239169.0)   ); // converted to J/mol from -57.16 kcal/mol as in table
        CHECK( props.H0  == Approx(-280617.0)   );
        CHECK( props.V0  == Approx(1.77559e-05) );
        CHECK( props.Cp0 == Approx(72.4980)     );
        CHECK( props.Cv0 == Approx(67.2659)     );
    }
}
