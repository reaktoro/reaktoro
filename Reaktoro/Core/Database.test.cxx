// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/Database.hpp>
using namespace Reaktoro;

namespace test {

auto createDatabase() -> Database
{
    Database db;

    // AQUEOUS SPECIES
    db.addSpecies( Species("H2O(aq)"  ) );
    db.addSpecies( Species("H+"       ) );
    db.addSpecies( Species("OH-"      ) );
    db.addSpecies( Species("H2(aq)"   ) );
    db.addSpecies( Species("O2(aq)"   ) );
    db.addSpecies( Species("Na+"      ) );
    db.addSpecies( Species("Cl-"      ) );
    db.addSpecies( Species("NaCl(aq)" ) );
    db.addSpecies( Species("HCl(aq)"  ) );
    db.addSpecies( Species("NaOH(aq)" ) );
    db.addSpecies( Species("Ca++"     ) );
    db.addSpecies( Species("Mg++"     ) );
    db.addSpecies( Species("CO2(aq)"  ) );
    db.addSpecies( Species("HCO3-"    ) );
    db.addSpecies( Species("CO3--"    ) );
    db.addSpecies( Species("CaCl2(aq)") );
    db.addSpecies( Species("MgCl2(aq)") );
    db.addSpecies( Species("SiO2(aq)" ) );

    // SOLID SPECIES
    db.addSpecies( Species("NaCl(s)"       ).withName("Halite")    );
    db.addSpecies( Species("CaCO3(s)"      ).withName("Calcite")   );
    db.addSpecies( Species("MgCO3(s)"      ).withName("Magnesite") );
    db.addSpecies( Species("CaMg(CO3)2(s)" ).withName("Dolomite")  );
    db.addSpecies( Species("SiO2(s)"       ).withName("Quartz")    );

    // GASEOUS SPECIES
    db.addSpecies( Species("CO2(g)" ) );
    db.addSpecies( Species("O2(g)"  ) );
    db.addSpecies( Species("H2(g)"  ) );
    db.addSpecies( Species("H2O(g)" ) );
    db.addSpecies( Species("CH4(g)" ) );
    db.addSpecies( Species("CO(g)"  ) );

    return db;
}

} // namespace test

TEST_CASE("Testing Database class", "[Database]")
{
    Database db = test::createDatabase();

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::elements
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.elements().index("H") );
    REQUIRE_NOTHROW( db.elements().index("C") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::species
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.species().index("H2O(aq)") );
    REQUIRE_NOTHROW( db.species().index("CO2(g)")  );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::speciesWithAggregateState
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H2O(aq)") );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H+")      );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Gas).index("CO2(g)")      );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Solid).index("Halite")    );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::addSpecies
    //-------------------------------------------------------------------------
    db.addSpecies( Species("Fe(s)").withName("Iron")   );
    db.addSpecies( Species("Cu(s)").withName("Copper") );

    REQUIRE_NOTHROW( db.species().index("Iron")   );
    REQUIRE_NOTHROW( db.species().index("Copper") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::attachData
    //-------------------------------------------------------------------------
    db.attachData( String("SomeData") );

    REQUIRE( std::any_cast<String>(db.attachedData()) == "SomeData" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::clear
    //-------------------------------------------------------------------------
    db.clear();

    REQUIRE( db.elements().size() == 0 );
    REQUIRE( db.species().size()  == 0 );
}
