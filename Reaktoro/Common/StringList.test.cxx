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
#include <Reaktoro/Common/StringList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StringList", "[StringList]")
{
    StringList list;

    // Test with constructor StringList(const std::char*)
    list = StringList("H2O H+ OH- Ca++ Ca(CO3)");

    REQUIRE(list.size() == 5);
    REQUIRE(list[0] == "H2O");
    REQUIRE(list[1] == "H+");
    REQUIRE(list[2] == "OH-");
    REQUIRE(list[3] == "Ca++");
    REQUIRE(list[4] == "Ca(CO3)");

    // Test with constructor StringList(const std::char*, const char token)
    list = StringList("H2O;H+;OH-;Ca++;Ca(CO3)", ';');

    REQUIRE(list.size() == 5);
    REQUIRE(list[0] == "H2O");
    REQUIRE(list[1] == "H+");
    REQUIRE(list[2] == "OH-");
    REQUIRE(list[3] == "Ca++");
    REQUIRE(list[4] == "Ca(CO3)");

    // Test with constructor StringList(const std::char*)
    list = StringList("hello there(spaces inside brackets are ignored!)");

    REQUIRE(list.size() == 2);
    REQUIRE(list[0] == "hello");
    REQUIRE(list[1] == "there(spaces inside brackets are ignored!)");

    // Test with constructor StringList(std::initializer_list<std::string>)
    list = StringList({ "H2O", "H+", "OH-", "Ca++", "Ca(CO3)" });

    REQUIRE(list.size() == 5);
    REQUIRE(list[0] == "H2O");
    REQUIRE(list[1] == "H+");
    REQUIRE(list[2] == "OH-");
    REQUIRE(list[3] == "Ca++");
    REQUIRE(list[4] == "Ca(CO3)");

    // Test with constructor StringList(std::vector<std::string>)
    std::vector<std::string> strs = { "H2O", "H+" };

    list = StringList(strs);

    REQUIRE(list.size() == 2);
    REQUIRE(list[0] == "H2O");
    REQUIRE(list[1] == "H+");
}
