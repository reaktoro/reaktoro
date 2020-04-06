// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
    Element element({
        .symbol = "H",
        .name = "Hydrogen",
        .atomic_number = 1,
        .atomic_weight = 0.001007940,
        .electronegativity = 2.20,
        .tags = {}
    });

    REQUIRE(element.symbol() == "H");
    REQUIRE(element.name() == "Hydrogen");
    REQUIRE(element.atomicNumber() == 1);
    REQUIRE(element.atomicWeight() == 0.001007940);
    REQUIRE(element.electronegativity() == 2.20);
    REQUIRE(element.molarMass() == element.atomicWeight());

    element = element.withSymbol("Na");
    REQUIRE(element.symbol() == "Na");

    element = element.withName("Sodium");
    REQUIRE(element.name() == "Sodium");

    element = element.withAtomicNumber(11);
    REQUIRE(element.atomicNumber() == 11);

    element = element.withAtomicWeight(0.022989768);
    REQUIRE(element.atomicWeight() == 0.022989768);

    element = element.withElectronegativity(0.93);
    REQUIRE(element.electronegativity() == 0.93);

    element = element.withTags({"tag1", "tag2"});
    REQUIRE(element.tags().size() == 2);
    REQUIRE(contains(element.tags(), "tag1"));
    REQUIRE(contains(element.tags(), "tag2"));
}
