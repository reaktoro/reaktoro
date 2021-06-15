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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/Element.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Element", "[Element]")
{
    Element element({ "H", 0.001007940, "Hydrogen" });

    REQUIRE(element.symbol() == "H");
    REQUIRE(element.name() == "Hydrogen");
    REQUIRE(element.molarMass() == 0.001007940);

    element = element.withSymbol("Na");
    REQUIRE(element.symbol() == "Na");

    element = element.withName("Sodium");
    REQUIRE(element.name() == "Sodium");

    element = element.withMolarMass(0.022989768);
    REQUIRE(element.molarMass() == 0.022989768);

    element = element.withTags({"tag1", "tag2"});
    REQUIRE(element.tags().size() == 2);
    REQUIRE(contains(element.tags(), "tag1"));
    REQUIRE(contains(element.tags(), "tag2"));
}
