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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelInterpolation class", "[StandardThermoModelInterpolation]")
{
    StandardThermoModelParamsInterpolation params;
    params.temperatures = {1.0, 2.0, 3.0};
    params.pressures = {4.0, 5.0};
    params.G0 = {{11.0, 12.0, 13.0},
                 {14.0, 15.0, 16.0}};
    params.H0 = {{21.0, 22.0, 23.0},
                 {24.0, 25.0, 26.0}};
    params.V0 = {{31.0, 32.0, 33.0},
                 {34.0, 35.0, 36.0}};
    params.Cp0 = {};
    params.VT0 = {};
    params.VP0 = {};
    params.Pref = 1.0;

    StandardThermoModel model = StandardThermoModelInterpolation(params);
    StandardThermoProps props;

    // CHECK WHEN (T, P) IS ON THE TOP-LEFT CORNER
    props = model(1.0, 4.0);

    CHECK( props.G0 == Approx(11.0 + props.V0*(4.0 - 1.0)) );
    CHECK( props.H0 == Approx(21.0) );
    CHECK( props.V0 == Approx(31.0) );
    CHECK( props.VT0 == 0.0 );
    CHECK( props.VP0 == 0.0 );
    CHECK( props.Cp0 == 0.0 );

    // CHECK WHEN (T, P) IS ON THE BOTTOM-RIGHT CORNER
    props = model(3.0, 5.0);

    CHECK( props.G0 == Approx(16.0 + props.V0*(5.0 - 1.0)) );
    CHECK( props.H0 == Approx(26.0) );
    CHECK( props.V0 == Approx(36.0) );
    CHECK( props.VT0 == 0.0 );
    CHECK( props.VP0 == 0.0 );
    CHECK( props.Cp0 == 0.0 );

    // CHECK WHEN (T, P) IS IN THE BOTTOM-RIGHT CELL (IN THE MIDDLE)
    props = model(2.5, 4.5);

    CHECK( props.G0 == Approx((12.0 + 13.0 + 15.0 + 16.0)/4 + props.V0*(4.5 - 1.0)) );
    CHECK( props.H0 == Approx((22.0 + 23.0 + 25.0 + 26.0)/4) );
    CHECK( props.V0 == Approx((32.0 + 33.0 + 35.0 + 36.0)/4) );
    CHECK( props.VT0 == 0.0 );
    CHECK( props.VP0 == 0.0 );
    CHECK( props.Cp0 == 0.0 );
}
