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
#include <Reaktoro/Core/PhaseList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing PhaseList", "[PhaseList]")
{
    PhaseList phases;
    PhaseList filtered;

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::append
    //-------------------------------------------------------------------------
    phases.append(Phase()
        .withName("AqueousPhase")
        .withSpecies(SpeciesList("H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq)"))
        .withStateOfMatter(StateOfMatter::Liquid));

    phases.append(Phase()
        .withName("GaseousPhase")
        .withSpecies(SpeciesList("H2O(g) CO2(g) CH4(g) O2(g) H2(g) CO(g)"))
        .withStateOfMatter(StateOfMatter::Gas));

    phases.append(Phase()
        .withName("Calcite")
        .withSpecies({ Species("CaCO3(s)").withName("Calcite") })
        .withStateOfMatter(StateOfMatter::Solid));

    phases.append(Phase()
        .withName("Halite")
        .withSpecies({ Species("NaCl(s)").withName("Halite") })
        .withStateOfMatter(StateOfMatter::Solid));

    phases.append(Phase()
        .withName("Magnesite")
        .withSpecies({ Species("MgCO3(s)").withName("Magnesite") })
        .withStateOfMatter(StateOfMatter::Solid));

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::size
    //-------------------------------------------------------------------------
    REQUIRE( phases.size() == 5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::operator[](Index)
    //-------------------------------------------------------------------------
    REQUIRE( phases[0].name() == "AqueousPhase" );
    REQUIRE( phases[1].name() == "GaseousPhase" );
    REQUIRE( phases[2].name() == "Calcite"         );
    REQUIRE( phases[3].name() == "Halite"          );
    REQUIRE( phases[4].name() == "Magnesite"       );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::find
    //-------------------------------------------------------------------------
    REQUIRE( phases.find("AqueousPhase") == 0 );
    REQUIRE( phases.find("GaseousPhase") == 1 );
    REQUIRE( phases.find("Calcite")         == 2 );
    REQUIRE( phases.find("Halite")          == 3 );
    REQUIRE( phases.find("Magnesite")       == 4 );

    REQUIRE( phases.find("@#$") >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( phases.findWithName("AqueousPhase") == 0 );
    REQUIRE( phases.findWithName("GaseousPhase") == 1 );
    REQUIRE( phases.findWithName("Calcite")         == 2 );
    REQUIRE( phases.findWithName("Halite")          == 3 );
    REQUIRE( phases.findWithName("Magnesite")       == 4 );

    REQUIRE( phases.findWithName("@#$") >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::findWithSpecies(index)
    //-------------------------------------------------------------------------
    REQUIRE( phases.findWithSpecies(0)  == 0 ); //  idx(0) -- H2O(aq)
    REQUIRE( phases.findWithSpecies(5)  == 0 ); //  idx(5) -- Na+
    REQUIRE( phases.findWithSpecies(8)  == 1 ); //  idx(8) -- H2O(g)
    REQUIRE( phases.findWithSpecies(13) == 1 ); // idx(13) -- CO(g)
    REQUIRE( phases.findWithSpecies(14) == 2 ); // idx(14) -- Calcite
    REQUIRE( phases.findWithSpecies(15) == 3 ); // idx(15) -- Halite
    REQUIRE( phases.findWithSpecies(16) == 4 ); // idx(16) -- Magnesite

    REQUIRE( phases.findWithSpecies(100) >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::findWithSpecies(name)
    //-------------------------------------------------------------------------
    REQUIRE( phases.findWithSpecies("H2O(aq)")   == 0 );
    REQUIRE( phases.findWithSpecies("Na+")       == 0 );
    REQUIRE( phases.findWithSpecies("H2O(g)")    == 1 );
    REQUIRE( phases.findWithSpecies("CO(g)")     == 1 );
    REQUIRE( phases.findWithSpecies("Calcite")   == 2 );
    REQUIRE( phases.findWithSpecies("Halite")    == 3 );
    REQUIRE( phases.findWithSpecies("Magnesite") == 4 );

    REQUIRE( phases.findWithSpecies("@#$") >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::findWithAggregateState
    //-------------------------------------------------------------------------
    REQUIRE( phases.findWithAggregateState(AggregateState::Aqueous) == 0 );
    REQUIRE( phases.findWithAggregateState(AggregateState::Gas)     == 1 );
    REQUIRE( phases.findWithAggregateState(AggregateState::Solid)   == 2 );

    REQUIRE( phases.findWithAggregateState(AggregateState::Undefined) >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::findWithStateOfMatter
    //-------------------------------------------------------------------------
    REQUIRE( phases.findWithStateOfMatter(StateOfMatter::Liquid) == 0 );
    REQUIRE( phases.findWithStateOfMatter(StateOfMatter::Gas)    == 1 );
    REQUIRE( phases.findWithStateOfMatter(StateOfMatter::Solid)  == 2 );

    REQUIRE( phases.findWithStateOfMatter(StateOfMatter::Plasma) >= phases.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::index
    //-------------------------------------------------------------------------
    REQUIRE( phases.index("AqueousPhase") == 0 );
    REQUIRE( phases.index("GaseousPhase") == 1 );
    REQUIRE( phases.index("Calcite")         == 2 );
    REQUIRE( phases.index("Halite")          == 3 );
    REQUIRE( phases.index("Magnesite")         == 4 );

    REQUIRE_THROWS( phases.index("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::indexWithName
    //-------------------------------------------------------------------------
    REQUIRE( phases.indexWithName("AqueousPhase") == 0 );
    REQUIRE( phases.indexWithName("GaseousPhase") == 1 );
    REQUIRE( phases.indexWithName("Calcite")         == 2 );
    REQUIRE( phases.indexWithName("Halite")          == 3 );
    REQUIRE( phases.indexWithName("Magnesite")         == 4 );

    REQUIRE_THROWS( phases.indexWithName("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::indexWithSpecies(index)
    //-------------------------------------------------------------------------
    REQUIRE( phases.indexWithSpecies(0)  == 0 ); //  idx(0) -- H2O(aq)
    REQUIRE( phases.indexWithSpecies(5)  == 0 ); //  idx(5) -- Na+
    REQUIRE( phases.indexWithSpecies(8)  == 1 ); //  idx(8) -- H2O(g)
    REQUIRE( phases.indexWithSpecies(13) == 1 ); // idx(13) -- CO(g)
    REQUIRE( phases.indexWithSpecies(14) == 2 ); // idx(14) -- Calcite
    REQUIRE( phases.indexWithSpecies(15) == 3 ); // idx(15) -- Halite
    REQUIRE( phases.indexWithSpecies(16) == 4 ); // idx(16) -- Magnesite

    REQUIRE_THROWS( phases.indexWithSpecies(100) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::indexWithSpecies(name)
    //-------------------------------------------------------------------------
    REQUIRE( phases.indexWithSpecies("H2O(aq)") == 0 );
    REQUIRE( phases.indexWithSpecies("Na+")     == 0 );
    REQUIRE( phases.indexWithSpecies("H2O(g)")  == 1 );
    REQUIRE( phases.indexWithSpecies("CO(g)")   == 1 );
    REQUIRE( phases.indexWithSpecies("Calcite") == 2 );
    REQUIRE( phases.indexWithSpecies("Halite")  == 3 );
    REQUIRE( phases.indexWithSpecies("Magnesite") == 4 );

    REQUIRE_THROWS( phases.indexWithSpecies("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::indexWithAggregateState
    //-------------------------------------------------------------------------
    REQUIRE( phases.indexWithAggregateState(AggregateState::Aqueous) == 0 );
    REQUIRE( phases.indexWithAggregateState(AggregateState::Gas)     == 1 );
    REQUIRE( phases.indexWithAggregateState(AggregateState::Solid)   == 2 );

    REQUIRE_THROWS( phases.indexWithAggregateState(AggregateState::Undefined) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::indexWithStateOfMatter
    //-------------------------------------------------------------------------
    REQUIRE( phases.indexWithStateOfMatter(StateOfMatter::Liquid) == 0 );
    REQUIRE( phases.indexWithStateOfMatter(StateOfMatter::Gas)    == 1 );
    REQUIRE( phases.indexWithStateOfMatter(StateOfMatter::Solid)  == 2 );

    REQUIRE_THROWS( phases.indexWithStateOfMatter(StateOfMatter::Plasma) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::withNames
    //-------------------------------------------------------------------------
    filtered = phases.withNames("AqueousPhase GaseousPhase");

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered[0].name() == "AqueousPhase" );
    REQUIRE( filtered[1].name() == "GaseousPhase" );

    REQUIRE_THROWS( filtered.withNames("Calcite @#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::withStateOfMatter
    //-------------------------------------------------------------------------
    filtered = phases.withStateOfMatter(StateOfMatter::Solid);

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "Calcite"   );
    REQUIRE( filtered[1].name() == "Halite"    );
    REQUIRE( filtered[2].name() == "Magnesite" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::withAggregateState
    //-------------------------------------------------------------------------
    filtered = phases.withAggregateState(AggregateState::Aqueous);

    REQUIRE( filtered.size() == 1 );
    REQUIRE( filtered[0].name() == "AqueousPhase" );

    filtered = phases.withAggregateState(AggregateState::Solid);

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "Calcite"   );
    REQUIRE( filtered[1].name() == "Halite"    );
    REQUIRE( filtered[2].name() == "Magnesite" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: PhaseList::numSpeciesUntilPhase
    //-------------------------------------------------------------------------
    REQUIRE( phases.numSpeciesUntilPhase(0) == 0 );
    REQUIRE( phases.numSpeciesUntilPhase(1) == 8 );
    REQUIRE( phases.numSpeciesUntilPhase(2) == 8 + 6 );
    REQUIRE( phases.numSpeciesUntilPhase(3) == 8 + 6 + 1 );
    REQUIRE( phases.numSpeciesUntilPhase(4) == 8 + 6 + 1 + 1 );
    REQUIRE( phases.numSpeciesUntilPhase(5) == 8 + 6 + 1 + 1 + 1 );
}
