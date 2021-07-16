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
#include <Reaktoro/Common/Units.hpp>

TEST_CASE("Testing Units module", "[Units]")
{
    auto check_identity = [](auto value, auto from, auto to)
    {
        REQUIRE( units::convert(convert(value, from, to), to, from) == Approx(value) );
    };

    auto x = GENERATE(0.0, 1.0, 100.0);

    INFO("x = " << x);

    //-------------------------------------------------------------------------
    // TEMPERATURE UNITS
    //-------------------------------------------------------------------------
    REQUIRE( units::convert(x, "degC", "K")    == Approx(x * 1.0 + 273.15) );
    REQUIRE( units::convert(x, "degC", "degF") == Approx(x * 1.8 +  32.00) );
    REQUIRE( units::convert(x, "degC", "degR") == Approx(x * 1.8 + 491.67) );

    REQUIRE( units::convert(x, "degC", "celsius")    == Approx(x).scale(1.0) );
    REQUIRE( units::convert(x, "degF", "fahrenheit") == Approx(x).scale(1.0) );
    REQUIRE( units::convert(x, "degR", "rankine")    == Approx(x).scale(1.0) );
    REQUIRE( units::convert(x, "K"   , "kelvin")     == Approx(x).scale(1.0) );

    //-------------------------------------------------------------------------
    // PRESSURES UNITS
    //-------------------------------------------------------------------------
    REQUIRE( units::convert(x, "Pa"     , "Pa") == Approx(x * 1.0)             );
    REQUIRE( units::convert(x, "kPa"    , "Pa") == Approx(x * 1.0e+3)          );
    REQUIRE( units::convert(x, "MPa"    , "Pa") == Approx(x * 1.0e+6)          );
    REQUIRE( units::convert(x, "GPa"    , "Pa") == Approx(x * 1.0e+9)          );
    REQUIRE( units::convert(x, "atm"    , "Pa") == Approx(x * 101325)          );
    REQUIRE( units::convert(x, "mmHg"   , "Pa") == Approx(x * 133.32239)       );
    REQUIRE( units::convert(x, "inHg"   , "Pa") == Approx(x * 3386.3886666667) );
    REQUIRE( units::convert(x, "psi"    , "Pa") == Approx(x * 6894.75729)      );
    REQUIRE( units::convert(x, "kpsi"   , "Pa") == Approx(x * 6894.75729e+3)   );
    REQUIRE( units::convert(x, "Mpsi"   , "Pa") == Approx(x * 6894.75729e+6)   );
    REQUIRE( units::convert(x, "psf"    , "Pa") == Approx(x * 47.88)           );
    REQUIRE( units::convert(x, "bar"    , "Pa") == Approx(x * 1.0e+5)          );
    REQUIRE( units::convert(x, "mbar"   , "Pa") == Approx(x * 1.0e+2)          );
    REQUIRE( units::convert(x, "kbar"   , "Pa") == Approx(x * 1.0e+8)          );
    REQUIRE( units::convert(x, "Mbar"   , "Pa") == Approx(x * 1.0e+11)         );
    REQUIRE( units::convert(x, "torr"   , "Pa") == Approx(x * 133.3223684211)  );
    REQUIRE( units::convert(x, "inH2O"  , "Pa") == Approx(x * 249.08891)       );
    REQUIRE( units::convert(x, "ftH2O"  , "Pa") == Approx(x * 249.08891 * 12)  );
    REQUIRE( units::convert(x, "pascal" , "Pa") == Approx(x * 1.0)             );
}
