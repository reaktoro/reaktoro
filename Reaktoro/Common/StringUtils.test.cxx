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
#include <Reaktoro/Common/StringUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StringUtils", "[StringUtils]")
{
    Strings words;

    REQUIRE( stringfy(".", "a", "b") == "a.b");
    REQUIRE( stringfy(".", 1, 2, "abc") == "1.2.abc");

    words = {"alpha", "beta", "gamma"};

    REQUIRE( str(words) == "alpha, beta, gamma");

    REQUIRE( str(1, 2, " abc") == "12 abc");

    REQUIRE( replace("CARBON DIOXIDE", " ", "-")   == "CARBON-DIOXIDE" );
    REQUIRE( replace("CARBON DIOXIDE", "O", "x")   == "CARBXN DIXXIDE" );
    REQUIRE( replace("CARBON DIOXIDE", "", "XYZ")  == "CARBON DIOXIDE" );

    words = {"alpha", "beta", "gamma", "beta", "gamma", "beta", "gamma", "alpha"};

    words = makeunique(words, "!");

    REQUIRE( words[0] == "alpha"   );
    REQUIRE( words[1] == "beta"    );
    REQUIRE( words[2] == "gamma"   );
    REQUIRE( words[3] == "beta!"   );
    REQUIRE( words[4] == "gamma!"  );
    REQUIRE( words[5] == "beta!!"  );
    REQUIRE( words[6] == "gamma!!" );
    REQUIRE( words[7] == "alpha!"  );
}
