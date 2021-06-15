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
#include <Reaktoro/Common/ParseUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ParseUtils module", "[ParseUtils]")
{
    //-----------------------------------------------------------------
    // Testing method: parseReaction
    //-----------------------------------------------------------------

    // TODO: Rename parseReaction to parseReactionEquation and implement tests for this method.

    //-----------------------------------------------------------------
    // Testing method: parseNumberStringPairs
    //-----------------------------------------------------------------
    Pairs<String, double> pairs;
    pairs = parseNumberStringPairs("2:H 1:O");
    CHECK( pairs == Pairs<String, double>{{"H", 2}, {"O", 1}} );

    pairs = parseNumberStringPairs("1:Ca 2:Cl");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"Cl", 2}} );

    pairs = parseNumberStringPairs("1:Mg 1:C 3:O");
    CHECK( pairs == Pairs<String, double>{{"Mg", 1}, {"C", 1}, {"O", 3}} );

    pairs = parseNumberStringPairs("1:Ca 1:Mg 2:C 6:O");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"Mg", 1}, {"C", 2}, {"O", 6}} );

    pairs = parseNumberStringPairs("3:Fe 2:Al 3:Si 12:O");
    CHECK( pairs == Pairs<String, double>{{"Fe", 3}, {"Al", 2}, {"Si", 3}, {"O", 12}} );

    pairs = parseNumberStringPairs("1:Fe 1:Al 1:Fe 3:Si 1:Fe 1:Al 4:O 4:O 4:O");
    CHECK( pairs == Pairs<String, double>{{"Fe", 3}, {"Al", 2}, {"Si", 3}, {"O", 12}} );

    pairs = parseNumberStringPairs("1:Na -1:E");
    CHECK( pairs == Pairs<String, double>{{"Na", 1}, {"E", -1}} );

    pairs = parseNumberStringPairs("1:Na -1:E -3:E 3:E");
    CHECK( pairs == Pairs<String, double>{{"Na", 1}, {"E", -1}} );

    pairs = parseNumberStringPairs("1:Ca -2:E");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"E", -2}} );

    pairs = parseNumberStringPairs("1:Fe -3:E");
    CHECK( pairs == Pairs<String, double>{{"Fe", 1}, {"E", -3}} );

    pairs = parseNumberStringPairs("1:C 3:O 2:E");
    CHECK( pairs == Pairs<String, double>{{"C", 1}, {"O", 3}, {"E", 2}} );

    pairs = parseNumberStringPairs("1:Ca++ 1:CO3--");
    CHECK( pairs == Pairs<String, double>{{"Ca++", 1}, {"CO3--", 1}} );

    pairs = parseNumberStringPairs("1:Ca++ 1:Mg++ 2:CO3--");
    CHECK( pairs == Pairs<String, double>{{"Ca++", 1}, {"Mg++", 1}, {"CO3--", 2}} );
}
