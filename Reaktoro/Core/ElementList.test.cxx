// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ElementList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ElementList", "[ElementList]")
{
    ElementList elements;
    ElementList filtered;

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ElementList(formulas)
    //-------------------------------------------------------------------------
    elements = {
        Element({ "H"  , "Hydrogen"      , 1   , 0.001007940 , 2.20, {"tag1"} }),
        Element({ "He" , "Helium"        , 2   , 0.004002602 , 0.00, {"tag1"} }),
        Element({ "Li" , "Lithium"       , 3   , 0.006941000 , 0.98, {"tag1"} }),
        Element({ "Be" , "Beryllium"     , 4   , 0.009012180 , 1.57, {"tag2"} }),
        Element({ "B"  , "Boron"         , 5   , 0.010811000 , 2.04, {"tag2"} }),
        Element({ "C"  , "Carbon"        , 6   , 0.012011000 , 2.55, {"tag2"} }),
        Element({ "N"  , "Nitrogen"      , 7   , 0.014006740 , 3.04, {"tag3"} }),
        Element({ "O"  , "Oxygen"        , 8   , 0.015999400 , 3.44, {"tag3"} }),
        Element({ "F"  , "Fluorine"      , 9   , 0.018998403 , 3.98, {"tag1", "tag3"} }),
    };

    REQUIRE( elements.size() == 9 );

    REQUIRE( elements[0].symbol()  == "H" );
    REQUIRE( elements[1].symbol()  == "He");
    REQUIRE( elements[2].symbol()  == "Li");
    REQUIRE( elements[3].symbol()  == "Be");
    REQUIRE( elements[4].symbol()  == "B" );
    REQUIRE( elements[5].symbol()  == "C" );
    REQUIRE( elements[6].symbol()  == "N" );
    REQUIRE( elements[7].symbol()  == "O" );
    REQUIRE( elements[8].symbol()  == "F" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::find
    //-------------------------------------------------------------------------
    REQUIRE( elements.find("H")  == 0);
    REQUIRE( elements.find("He") == 1);
    REQUIRE( elements.find("Li") == 2);
    REQUIRE( elements.find("Be") == 3);
    REQUIRE( elements.find("B")  == 4);
    REQUIRE( elements.find("C")  == 5);
    REQUIRE( elements.find("N")  == 6);
    REQUIRE( elements.find("O")  == 7);
    REQUIRE( elements.find("F")  == 8);

    REQUIRE( elements.find("Xy") >= elements.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::findWithSymbol
    //-------------------------------------------------------------------------
    REQUIRE( elements.findWithSymbol("H")  == 0);
    REQUIRE( elements.findWithSymbol("He") == 1);
    REQUIRE( elements.findWithSymbol("Li") == 2);
    REQUIRE( elements.findWithSymbol("Be") == 3);
    REQUIRE( elements.findWithSymbol("B")  == 4);
    REQUIRE( elements.findWithSymbol("C")  == 5);
    REQUIRE( elements.findWithSymbol("N")  == 6);
    REQUIRE( elements.findWithSymbol("O")  == 7);
    REQUIRE( elements.findWithSymbol("F")  == 8);

    REQUIRE( elements.findWithSymbol("Xy") >= elements.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( elements.findWithName("Hydrogen")  == 0);
    REQUIRE( elements.findWithName("Helium")    == 1);
    REQUIRE( elements.findWithName("Lithium")   == 2);
    REQUIRE( elements.findWithName("Beryllium") == 3);
    REQUIRE( elements.findWithName("Boron")     == 4);
    REQUIRE( elements.findWithName("Carbon")    == 5);
    REQUIRE( elements.findWithName("Nitrogen")  == 6);
    REQUIRE( elements.findWithName("Oxygen")    == 7);
    REQUIRE( elements.findWithName("Fluorine")  == 8);

    REQUIRE( elements.findWithName("Xyrium") >= elements.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::index
    //-------------------------------------------------------------------------
    REQUIRE( elements.index("H")  == 0);
    REQUIRE( elements.index("He") == 1);
    REQUIRE( elements.index("Li") == 2);
    REQUIRE( elements.index("Be") == 3);
    REQUIRE( elements.index("B")  == 4);
    REQUIRE( elements.index("C")  == 5);
    REQUIRE( elements.index("N")  == 6);
    REQUIRE( elements.index("O")  == 7);
    REQUIRE( elements.index("F")  == 8);

    REQUIRE_THROWS( elements.index("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::indexWithSymbol
    //-------------------------------------------------------------------------
    REQUIRE( elements.indexWithSymbol("H")  == 0);
    REQUIRE( elements.indexWithSymbol("He") == 1);
    REQUIRE( elements.indexWithSymbol("Li") == 2);
    REQUIRE( elements.indexWithSymbol("Be") == 3);
    REQUIRE( elements.indexWithSymbol("B")  == 4);
    REQUIRE( elements.indexWithSymbol("C")  == 5);
    REQUIRE( elements.indexWithSymbol("N")  == 6);
    REQUIRE( elements.indexWithSymbol("O")  == 7);
    REQUIRE( elements.indexWithSymbol("F")  == 8);

    REQUIRE_THROWS( elements.indexWithSymbol("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::indexWithName
    //-------------------------------------------------------------------------
    REQUIRE( elements.indexWithName("Hydrogen")  == 0);
    REQUIRE( elements.indexWithName("Helium")    == 1);
    REQUIRE( elements.indexWithName("Lithium")   == 2);
    REQUIRE( elements.indexWithName("Beryllium") == 3);
    REQUIRE( elements.indexWithName("Boron")     == 4);
    REQUIRE( elements.indexWithName("Carbon")    == 5);
    REQUIRE( elements.indexWithName("Nitrogen")  == 6);
    REQUIRE( elements.indexWithName("Oxygen")    == 7);
    REQUIRE( elements.indexWithName("Fluorine")  == 8);

    REQUIRE_THROWS( elements.indexWithName("Xyrium") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::get
    //-------------------------------------------------------------------------
    REQUIRE( elements.get("H").symbol()  == "H"  );
    REQUIRE( elements.get("He").symbol() == "He" );
    REQUIRE( elements.get("Li").symbol() == "Li" );
    REQUIRE( elements.get("Be").symbol() == "Be" );
    REQUIRE( elements.get("B").symbol()  == "B"  );
    REQUIRE( elements.get("C").symbol()  == "C"  );
    REQUIRE( elements.get("N").symbol()  == "N"  );
    REQUIRE( elements.get("O").symbol()  == "O"  );
    REQUIRE( elements.get("F").symbol()  == "F"  );

    REQUIRE_THROWS( elements.get("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::getWithSymbol
    //-------------------------------------------------------------------------
    REQUIRE( elements.getWithSymbol("H") .symbol() == "H"  );
    REQUIRE( elements.getWithSymbol("He").symbol() == "He" );
    REQUIRE( elements.getWithSymbol("Li").symbol() == "Li" );
    REQUIRE( elements.getWithSymbol("Be").symbol() == "Be" );
    REQUIRE( elements.getWithSymbol("B") .symbol() == "B"  );
    REQUIRE( elements.getWithSymbol("C") .symbol() == "C"  );
    REQUIRE( elements.getWithSymbol("N") .symbol() == "N"  );
    REQUIRE( elements.getWithSymbol("O") .symbol() == "O"  );
    REQUIRE( elements.getWithSymbol("F") .symbol() == "F"  );

    REQUIRE_THROWS( elements.getWithSymbol("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::getWithName
    //-------------------------------------------------------------------------
    REQUIRE( elements.getWithName("Hydrogen") .name() == "Hydrogen" );
    REQUIRE( elements.getWithName("Helium")   .name() == "Helium"   );
    REQUIRE( elements.getWithName("Lithium")  .name() == "Lithium"  );
    REQUIRE( elements.getWithName("Beryllium").name() == "Beryllium");
    REQUIRE( elements.getWithName("Boron")    .name() == "Boron"    );
    REQUIRE( elements.getWithName("Carbon")   .name() == "Carbon"   );
    REQUIRE( elements.getWithName("Nitrogen") .name() == "Nitrogen" );
    REQUIRE( elements.getWithName("Oxygen")   .name() == "Oxygen"   );
    REQUIRE( elements.getWithName("Fluorine") .name() == "Fluorine" );

    REQUIRE_THROWS( elements.getWithName("Xyrium") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withSymbols
    //-------------------------------------------------------------------------
    filtered = elements.withSymbols("He O C");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered.indexWithSymbol("He") < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("O" ) < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("C" ) < filtered.size() );

    REQUIRE_THROWS( elements.withSymbols("He ABC O XYZ C") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withNames
    //-------------------------------------------------------------------------
    filtered = elements.withNames("Helium Oxygen Carbon");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered.indexWithName("Helium") < filtered.size() );
    REQUIRE( filtered.indexWithName("Oxygen") < filtered.size() );
    REQUIRE( filtered.indexWithName("Carbon") < filtered.size() );

    REQUIRE_THROWS( elements.withSymbols("Helium Xylium") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withTag
    //-------------------------------------------------------------------------
    filtered = elements.withTag("tag1");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithSymbol("H" ) < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("He") < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("Li") < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("F" ) < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withoutTag
    //-------------------------------------------------------------------------
    filtered = elements.withoutTag("tag1");

    REQUIRE( filtered.size() == 5 );
    REQUIRE( filtered.indexWithSymbol("Be") < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("B" ) < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("C" ) < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("N" ) < filtered.size() );
    REQUIRE( filtered.indexWithSymbol("O" ) < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withTags
    //-------------------------------------------------------------------------
    filtered = elements.withTags({"tag1", "tag2"});
    REQUIRE( filtered.size() == 0 );

    filtered = elements.withTags({"tag1", "tag3"});
    REQUIRE( filtered.size() == 1 );
    REQUIRE( filtered.index("F") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::withoutTags
    //-------------------------------------------------------------------------
    filtered = elements.withTags({"tag1", "tag2", "tag3"});

    REQUIRE( filtered.size() == 0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ElementList::append
    //-------------------------------------------------------------------------
    elements.append(Element().withSymbol("Xy").withName("Xyrium"));

    REQUIRE( elements.indexWithSymbol("Xy")   < elements.size() );
    REQUIRE( elements.indexWithName("Xyrium") < elements.size() );
}
