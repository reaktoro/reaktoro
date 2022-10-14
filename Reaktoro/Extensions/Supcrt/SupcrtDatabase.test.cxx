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
#include <Reaktoro/Extensions/Supcrt/SupcrtDatabase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SupcrtDatabase object creation using embedded files", "[SupcrtDatabase]")
{
    CHECK_NOTHROW( SupcrtDatabase("supcrt98") );
    CHECK_NOTHROW( SupcrtDatabase("supcrt98-organics") );
    CHECK_NOTHROW( SupcrtDatabase("supcrt07") );
    CHECK_NOTHROW( SupcrtDatabase("supcrt07-organics") );
    CHECK_NOTHROW( SupcrtDatabase("supcrt16") );
    CHECK_NOTHROW( SupcrtDatabase("supcrt16-organics") );
    CHECK_NOTHROW( SupcrtDatabase("supcrtbl") );
    CHECK_NOTHROW( SupcrtDatabase("supcrtbl-organics") );
}

TEST_CASE("Testing SupcrtDatabase module", "[SupcrtDatabase]")
{
    SupcrtDatabase db("supcrt98");

    CHECK( db.species().findWithName("H2O(aq)") );
    CHECK( db.species().findWithName("H+") );
    CHECK( db.species().findWithName("OH-") );
    CHECK( db.species().findWithName("CO2(aq)") );

    // TODO: More tests are needed for SupcrtDatabase, at different T and P values,
    {
        auto sp = db.species().getWithName("H2O(aq)");

        const auto T = 25.0 + 273.15;
        const auto P =  1.0 * 1e5;

        auto [G0, H0, V0, Cp0, VT0, VP0] = sp.standardThermoProps(T, P);

        CHECK( G0  == Approx(-237182)     );
        CHECK( H0  == Approx(-285831)     );
        CHECK( V0  == Approx(1.80686e-05) );
        CHECK( Cp0 == Approx(75.3276)     );
    }

    {
        auto sp = db.species().getWithName("HCO3-");

        const auto T = 25.0 + 273.15;
        const auto P =  1.0 * 1e5;

        auto [G0, H0, V0, Cp0, VT0, VP0] = sp.standardThermoProps(T, P);

        CHECK( G0  == Approx(-586940)     );
        CHECK( H0  == Approx(-689934)     );
        CHECK( V0  == Approx(2.42086e-05) );
        CHECK( Cp0 == Approx(-34.9279)    );
    }

    {
        auto sp = db.species().getWithName("CO2(g)");

        const auto T = 25.0 + 273.15;
        const auto P =  1.0 * 1e5;

        auto [G0, H0, V0, Cp0, VT0, VP0] = sp.standardThermoProps(T, P);

        CHECK( G0  == Approx(-394359)     );
        CHECK( H0  == Approx(-393509)     );
        CHECK( V0  == Approx(0.0).scale(1.0)   );
        CHECK( Cp0 == Approx(37.1486)     );
    }

    {
        auto sp = db.species().getWithName("Calcite");

        const auto T = 25.0 + 273.15;
        const auto P =  1.0 * 1e5;

        auto [G0, H0, V0, Cp0, VT0, VP0] = sp.standardThermoProps(T, P);

        CHECK( G0  == Approx(-1.12918e+06));
        CHECK( H0  == Approx(-1.2073e+06) );
        CHECK( V0  == Approx(3.6934e-05)  );
        CHECK( Cp0 == Approx(81.8711)     );
    }
}
