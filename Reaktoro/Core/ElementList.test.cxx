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
#include <Reaktoro/Core/ElementList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ElementList", "[ElementList]")
{
    ElementList elements;
    ElementList filtered;

    // Test constructor ElementList(formulas)
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

    // Test method ElementList::indexWithSymbol
    REQUIRE( elements.indexWithSymbol("H")  == 0);
    REQUIRE( elements.indexWithSymbol("He") == 1);
    REQUIRE( elements.indexWithSymbol("Li") == 2);
    REQUIRE( elements.indexWithSymbol("Be") == 3);
    REQUIRE( elements.indexWithSymbol("B")  == 4);
    REQUIRE( elements.indexWithSymbol("C")  == 5);
    REQUIRE( elements.indexWithSymbol("N")  == 6);
    REQUIRE( elements.indexWithSymbol("O")  == 7);
    REQUIRE( elements.indexWithSymbol("F")  == 8);

    // Test method ElementList::indexWithName
    REQUIRE( elements.indexWithName("Hydrogen")  == 0);
    REQUIRE( elements.indexWithName("Helium")    == 1);
    REQUIRE( elements.indexWithName("Lithium")   == 2);
    REQUIRE( elements.indexWithName("Beryllium") == 3);
    REQUIRE( elements.indexWithName("Boron")     == 4);
    REQUIRE( elements.indexWithName("Carbon")    == 5);
    REQUIRE( elements.indexWithName("Nitrogen")  == 6);
    REQUIRE( elements.indexWithName("Oxygen")    == 7);
    REQUIRE( elements.indexWithName("Fluorine")  == 8);

    // Test method ElementList::withSymbols
    filtered = elements.withSymbols("He ABC O XYZ C");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].symbol() == "He" );
    REQUIRE( filtered[1].symbol() == "O"  );
    REQUIRE( filtered[2].symbol() == "C"  );

    // Test method ElementList::withNames
    filtered = elements.withNames("Helium ABC Oxygen XYZ Carbon");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "Helium" );
    REQUIRE( filtered[1].name() == "Oxygen" );
    REQUIRE( filtered[2].name() == "Carbon" );

    // Test method ElementList::withTag
    filtered = elements.withTag("tag1");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered[0].symbol()  == "H"  );
    REQUIRE( filtered[1].symbol()  == "He" );
    REQUIRE( filtered[2].symbol()  == "Li" );
    REQUIRE( filtered[3].symbol()  == "F"  );

    // Test method ElementList::withoutTag
    filtered = elements.withoutTag("tag1");

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered[0].symbol()  == "H"  );
    REQUIRE( filtered[1].symbol()  == "He" );
    REQUIRE( filtered[2].symbol()  == "Li" );
    REQUIRE( filtered[3].symbol()  == "Be" );
    REQUIRE( filtered[4].symbol()  == "B"  );
    REQUIRE( filtered[5].symbol()  == "C"  );

    // Test method ElementList::withTags
    filtered = elements.withTags({"tag1", "tag2"});

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered[0].symbol()  == "H"  );
    REQUIRE( filtered[1].symbol()  == "He" );
    REQUIRE( filtered[2].symbol()  == "Li" );
    REQUIRE( filtered[3].symbol()  == "Be" );
    REQUIRE( filtered[4].symbol()  == "B"  );
    REQUIRE( filtered[5].symbol()  == "C"  );

    // Test method ElementList::withoutTags
    filtered = elements.withTags({"tag1", "tag2", "tag3"});

    REQUIRE( filtered.size() == 0 );

    // Test method ElementList::append
    elements.append(Element("Xy").withName("Xyrium"));

    REQUIRE( elements.indexWithSymbol("Xy")   < elements.size() );
    REQUIRE( elements.indexWithName("Xyrium") < elements.size() );
}
