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
#include <Reaktoro/Water/WaterInterpolation.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>
#include <Reaktoro/Water/WaterUtils.hpp>
using namespace Reaktoro;

#define CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP(T, PMPa, reltol)                          \
    {                                                                                     \
        auto expected = waterDensityWagnerPruss(T, PMPa*1e6).val();                       \
        auto actual = waterDensityWagnerPrussInterp(T, PMPa*1e6).val();                   \
        CHECK( actual == Approx(expected).epsilon(reltol) );                              \
    }

#define CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP(T, PMPa, reltol)         \
    {                                                                         \
        auto expected = waterThermoPropsWagnerPruss(T, PMPa*1e6).D.val();     \
        auto actual = waterThermoPropsWagnerPrussInterp(T, PMPa*1e6).D.val(); \
        CHECK( actual == Approx(expected).epsilon(reltol) );                  \
    }

TEST_CASE("Testing water interpolation methods", "[WaterInterpolation]")
{
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 300.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 354.467,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 354.468,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 400.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP(1273.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP(1273.000, 1000.00, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 353.000,    2.00, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 423.000,    8.00, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 566.000,   14.00, 0.01); // permit 1% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 690.000,   87.00, 0.02); // permit 2% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 723.000,  125.00, 0.05); // permit 5% error deviation
    CHECK_WATER_DENSITY_WAGNER_PRUSS_INTERP( 834.000,  345.00, 0.02); // permit 2% error deviation

    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 300.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 354.467,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 354.468,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 400.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP(1273.000,    0.05, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP(1273.000, 1000.00, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 353.000,    2.00, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 423.000,    8.00, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 566.000,   14.00, 0.01); // permit 1% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 690.000,   87.00, 0.02); // permit 2% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 723.000,  125.00, 0.05); // permit 5% error deviation
    CHECK_WATER_THERMO_PROPS_WAGNER_PRUSS_INTERP( 834.000,  345.00, 0.02); // permit 2% error deviation

    // TODO: To reduce errors, above (note the 4.17% error at 723K and 125MPa), more refinement in the grid is needed.
}
