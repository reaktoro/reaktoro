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
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/ReactionList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionList", "[ReactionList]")
{
    SpeciesList species("H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) NaCl(s) O2(g)");
    Database db(species);

    ReactionList reactions;
    ReactionList filtered;

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::append
    //-------------------------------------------------------------------------
    reactions.append(db.reaction("H2O(aq) = H+ + OH-"));
    reactions.append(db.reaction("H2O(aq) = H2(aq) + 0.5*O2(aq)"));
    reactions.append(db.reaction("O2(g) = O2(aq)"));
    reactions.append(db.reaction("NaCl(s) = Na+ + Cl-"));
    reactions.append(db.reaction("NaCl(s)"));

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::size
    //-------------------------------------------------------------------------
    REQUIRE( reactions.size() == 5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::operator[](Index)
    //-------------------------------------------------------------------------
    REQUIRE( reactions[0].name() == "H2O(aq) = H+ + OH-" );
    REQUIRE( reactions[1].name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );
    REQUIRE( reactions[2].name() == "O2(g) = O2(aq)" );
    REQUIRE( reactions[3].name() == "NaCl(s) = Na+ + Cl-" );
    REQUIRE( reactions[4].name() == "NaCl(s)" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::find
    //-------------------------------------------------------------------------
    REQUIRE( reactions.find("H2O(aq) = H+ + OH-")            == 0 );
    REQUIRE( reactions.find("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1 );
    REQUIRE( reactions.find("O2(g) = O2(aq)")                == 2 );
    REQUIRE( reactions.find("NaCl(s) = Na+ + Cl-")           == 3 );
    REQUIRE( reactions.find("NaCl(s)")                       == 4 );

    REQUIRE( reactions.find("@#$") >= reactions.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( reactions.findWithName("H2O(aq) = H+ + OH-")            == 0 );
    REQUIRE( reactions.findWithName("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1 );
    REQUIRE( reactions.findWithName("O2(g) = O2(aq)")                == 2 );
    REQUIRE( reactions.findWithName("NaCl(s) = Na+ + Cl-")           == 3 );
    REQUIRE( reactions.findWithName("NaCl(s)")                       == 4 );

    REQUIRE( reactions.findWithName("@#$") >= reactions.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::index
    //-------------------------------------------------------------------------
    REQUIRE( reactions.index("H2O(aq) = H+ + OH-")            == 0 );
    REQUIRE( reactions.index("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1 );
    REQUIRE( reactions.index("O2(g) = O2(aq)")                == 2 );
    REQUIRE( reactions.index("NaCl(s) = Na+ + Cl-")           == 3 );
    REQUIRE( reactions.index("NaCl(s)")                       == 4 );

    REQUIRE_THROWS( reactions.index("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::indexWithName
    //-------------------------------------------------------------------------
    REQUIRE( reactions.indexWithName("H2O(aq) = H+ + OH-")            == 0 );
    REQUIRE( reactions.indexWithName("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1 );
    REQUIRE( reactions.indexWithName("O2(g) = O2(aq)")                == 2 );
    REQUIRE( reactions.indexWithName("NaCl(s) = Na+ + Cl-")           == 3 );
    REQUIRE( reactions.indexWithName("NaCl(s)")                       == 4 );

    REQUIRE_THROWS( reactions.indexWithName("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::get
    //-------------------------------------------------------------------------
    REQUIRE( reactions.get("H2O(aq) = H+ + OH-").name()            == "H2O(aq) = H+ + OH-" );
    REQUIRE( reactions.get("H2O(aq) = H2(aq) + 0.5*O2(aq)").name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );
    REQUIRE( reactions.get("O2(g) = O2(aq)").name()                == "O2(g) = O2(aq)" );
    REQUIRE( reactions.get("NaCl(s) = Na+ + Cl-").name()           == "NaCl(s) = Na+ + Cl-" );
    REQUIRE( reactions.get("NaCl(s)").name()                       == "NaCl(s)" );

    REQUIRE_THROWS( reactions.get("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::getWithName
    //-------------------------------------------------------------------------
    REQUIRE( reactions.getWithName("H2O(aq) = H+ + OH-").name()            == "H2O(aq) = H+ + OH-" );
    REQUIRE( reactions.getWithName("H2O(aq) = H2(aq) + 0.5*O2(aq)").name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );
    REQUIRE( reactions.getWithName("O2(g) = O2(aq)").name()                == "O2(g) = O2(aq)" );
    REQUIRE( reactions.getWithName("NaCl(s) = Na+ + Cl-").name()           == "NaCl(s) = Na+ + Cl-" );
    REQUIRE( reactions.getWithName("NaCl(s)").name()                       == "NaCl(s)" );

    REQUIRE_THROWS( reactions.getWithName("@#$") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ReactionList::withNames
    //-------------------------------------------------------------------------
    filtered = reactions.withNames({"H2O(aq) = H+ + OH-", "H2O(aq) = H2(aq) + 0.5*O2(aq)"});

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered[0].name() == "H2O(aq) = H+ + OH-" );
    REQUIRE( filtered[1].name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );

    REQUIRE_THROWS( filtered.withNames("Calcite @#$") );
}
