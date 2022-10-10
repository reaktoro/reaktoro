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
#include <Reaktoro/Water/WaterThermoProps.hpp>
using namespace Reaktoro;

#define CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, value) \
    CHECK( props.T.val()   == Approx(value) );            \
    CHECK( props.V.val()   == Approx(value) );            \
    CHECK( props.S.val()   == Approx(value) );            \
    CHECK( props.A.val()   == Approx(value) );            \
    CHECK( props.U.val()   == Approx(value) );            \
    CHECK( props.H.val()   == Approx(value) );            \
    CHECK( props.G.val()   == Approx(value) );            \
    CHECK( props.Cv.val()  == Approx(value) );            \
    CHECK( props.Cp.val()  == Approx(value) );            \
    CHECK( props.D.val()   == Approx(value) );            \
    CHECK( props.DT.val()  == Approx(value) );            \
    CHECK( props.DP.val()  == Approx(value) );            \
    CHECK( props.DTT.val() == Approx(value) );            \
    CHECK( props.DTP.val() == Approx(value) );            \
    CHECK( props.DPP.val() == Approx(value) );            \
    CHECK( props.P.val()   == Approx(value) );            \
    CHECK( props.PT.val()  == Approx(value) );            \
    CHECK( props.PD.val()  == Approx(value) );            \
    CHECK( props.PTT.val() == Approx(value) );            \
    CHECK( props.PTD.val() == Approx(value) );            \
    CHECK( props.PDD.val() == Approx(value) );


TEST_CASE("Testing class WaterThermoProps", "[WaterThermoProps]")
{
    const WaterThermoProps a = {
        1.0,  // T
        2.0,  // V
        3.0,  // S
        4.0,  // A
        5.0,  // U
        6.0,  // H
        7.0,  // G
        8.0,  // Cv
        9.0,  // Cp
        10.0, // D
        11.0, // DT
        12.0, // DP
        13.0, // DTT
        14.0, // DTP
        15.0, // DPP
        16.0, // P
        17.0, // PT
        18.0, // PD
        19.0, // PTT
        20.0, // PTD
        21.0, // PDD
    };

    const WaterThermoProps b = {
        -1.0,  // T
        -2.0,  // V
        -3.0,  // S
        -4.0,  // A
        -5.0,  // U
        -6.0,  // H
        -7.0,  // G
        -8.0,  // Cv
        -9.0,  // Cp
        -10.0, // D
        -11.0, // DT
        -12.0, // DP
        -13.0, // DTT
        -14.0, // DTP
        -15.0, // DPP
        -16.0, // P
        -17.0, // PT
        -18.0, // PD
        -19.0, // PTT
        -20.0, // PTD
        -21.0, // PDD
    };

    const WaterThermoProps c = {
        10.0, // T
        10.0, // V
        10.0, // S
        10.0, // A
        10.0, // U
        10.0, // H
        10.0, // G
        10.0, // Cv
        10.0, // Cp
        10.0, // D
        10.0, // DT
        10.0, // DP
        10.0, // DTT
        10.0, // DTP
        10.0, // DPP
        10.0, // P
        10.0, // PT
        10.0, // PD
        10.0, // PTT
        10.0, // PTD
        10.0, // PDD
    };

    WaterThermoProps props;

    props = a;
    props += b;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 0.0);

    props = a;
    props -= a;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 0.0);

    props = c;
    props *= 2.0;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 20.0);

    props = c;
    props *= real(3.0);
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 30.0);

    props = c;
    props /= 2.0;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 5.0);

    props = c;
    props /= real(5.0);
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 2.0);

    props = a + b;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 0.0);

    props = (-b) - (+a);
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 0.0);

    props = 2.0 * c;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 20.0);

    props = c * 3.0;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 30.0);

    props = 4.0 * c;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 40.0);

    props = c * 5.0;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 50.0);

    props = real(4.0) * c;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 40.0);

    props = c * real(5.0);
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 50.0);

    props = c / 5.0;
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 2.0);

    props = c / real(2.0);
    CHECK_WATER_THERMO_PROPS_VALUES_ARE(props, 5.0);
}
