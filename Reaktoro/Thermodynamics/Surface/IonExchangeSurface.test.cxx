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
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>
#include <Reaktoro/Singletons/Elements.hpp>

using namespace Reaktoro;

namespace test { extern auto createDatabase() -> Database; }

TEST_CASE("Testing IonExchangeSurface", "[IonExchangeSurface]")
{
    // Extend the list of elements by the exchanger element
    Elements::append(Element().withSymbol("X").withMolarMass(10.0));

    // Create custom database
    Database db = test::createDatabase();

    // Define ion exchange species list
    // Expected species: X- AlX3 CaX2 KX MgX2 NaX NH4X
    SpeciesList species_db = db.species().withAggregateState(AggregateState::IonExchange);

    // Create the aqueous mixture
    IonExchangeSurface surface_db(species_db);

    SECTION("Checking the charges of the species in IonExchangeSurface (custom database)")
    {
        // The numbers of exchanger's equivalents for exchange species
        ArrayXd ze = surface_db.ze();

        CHECK( ze[0] == 0 ); // X-
        CHECK( ze[1] == 3 ); // AlX3
        CHECK( ze[2] == 2 ); // CaX2
        CHECK( ze[3] == 1 ); // KX
        CHECK( ze[4] == 2 ); // MgX2
        CHECK( ze[5] == 1 ); // NaX
        CHECK( ze[6] == 1 ); // NH4X
    }

    SECTION("Checking indices in IonExchangeSurface (custom database)")
    {
        auto exchanger_index  = surface_db.indexExchanger();
        auto exchange_indices = surface_db.indicesExchange();

        CHECK( exchanger_index     == 0 );
        CHECK( exchange_indices[0] == 1 );
        CHECK( exchange_indices[1] == 2 );
        CHECK( exchange_indices[2] == 3 );
        CHECK( exchange_indices[3] == 4 );
        CHECK( exchange_indices[4] == 5 );
        CHECK( exchange_indices[5] == 6 );
    }

    SECTION("Checking the species in IonExchangeSurface (custom database)")
    {
        CHECK(surface_db.species()[0].name() == "X-"   ); // X-
        CHECK(surface_db.species()[1].name() == "AlX3" ); // AlX3
        CHECK(surface_db.species()[2].name() == "CaX2" ); // CaX2
        CHECK(surface_db.species()[3].name() == "KX"   ); // KX
        CHECK(surface_db.species()[4].name() == "MgX2" ); // MgX2
        CHECK(surface_db.species()[5].name() == "NaX"  ); // NaX
        CHECK(surface_db.species()[6].name() == "NH4X" ); // NH4X
    }
}
