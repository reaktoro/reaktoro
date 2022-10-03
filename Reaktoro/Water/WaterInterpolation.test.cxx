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
#include <Reaktoro/Water/WaterUtils.hpp>
using namespace Reaktoro;

#define CHECK_DENSITY_LIQUID(T, PMPa) \
    CHECK( waterDensityWagnerPrussInterp(T, PMPa*1e6) == Approx(waterDensityWagnerPruss(T, PMPa*1e6, StateOfMatter::Liquid)) );

#define CHECK_DENSITY_GAS(T, PMPa) \
    CHECK( waterDensityWagnerPrussInterp(T, PMPa*1e6) == Approx(waterDensityWagnerPruss(T, PMPa*1e6, StateOfMatter::Gas)) );

TEST_CASE("Testing water interpolation methods", "[WaterInterpolation]")
{
    const auto MPa = 1e6;

    CHECK( waterDensityWagnerPrussInterp( 300.000, 0.05*MPa).val()     == Approx(996.534) );      // error 1.093747046920056e-05 % when compared to 996.5338910044 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp( 354.467, 0.05*MPa).val()     == Approx(970.942) );      // error 1.017359326934713e-05 % when compared to 970.9420987797 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp( 354.468, 0.05*MPa).val()     == Approx(0.3086359765) ); // error 0.000774660655486% when compared to 0.3086383674 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp( 400.000, 0.05*MPa).val()     == Approx(0.27229) );      // error 0.001559709779697% when compared to 0.272294247 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(1273.000, 0.05*MPa).val()     == Approx(0.08511) );      // error 0.003205366516750% when compared to 0.085107272 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(1273.000, 1000*MPa).val()     == Approx(809.28) );       // error 3.282837988918427e-05 % when compared to 809.2802656736 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(353.000,   2*MPa).val()       == Approx(972.73352) );    // error 5.064858354630614e-05 % when compared to 972.7330273245 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(423.000,   8*MPa).val()       == Approx(921.35355) );    // error 0.000772530481138% when compared to 921.360667792 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(566.000,  14*MPa).val()       == Approx(738.0114944) );  // error 0.022461362711453% when compared to 738.1772990806 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(690.000,  87*MPa).val()       == Approx(650.7964578) );  // error 1.072531461874359% when compared to 643.8905292933 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(723.000, 125*MPa).val()       == Approx(686.10870775) ); // error 4.170730301765396% when compared to 658.6386653549 from waterDensityWagnerPruss
    CHECK( waterDensityWagnerPrussInterp(834.000, 345*MPa).val()       == Approx(757.17433517) ); // error 0.520612296527317% when compared to 753.2528084254 from waterDensityWagnerPruss

    CHECK( waterDensityWagnerPrussInterp( 370.000, 0.15*MPa).val()     == Approx(waterDensityWagnerPruss( 370.000, 0.15*MPa, StateOfMatter::Liquid)) );

    // TODO: To reduce errors, above (note the 4.17% error at 723K and 125MPa), more refinement in the grid is needed.
}
