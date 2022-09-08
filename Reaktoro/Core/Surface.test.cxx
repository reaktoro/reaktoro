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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/Surface.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Surface", "[Surface]")
{
    const auto phase1 = Phase()
        .withName("AqueousPhase")
        .withSpecies(SpeciesList("H2O(aq) H+(aq) OH-(aq) CO2(aq) HCO3-(aq) CO3--(aq)"))
        .withStateOfMatter(StateOfMatter::Liquid)
        ;

    const auto phase2 = Phase()
        .withName("GaseousPhase")
        .withSpecies(SpeciesList("CO2(g) H2O(g)"))
        .withStateOfMatter(StateOfMatter::Gas)
        ;

    WHEN("using constructor Surface()")
    {
        Surface surface;

        surface = surface.withName("AqueousPhase:GaseousPhase");
        CHECK(surface.name() == "AqueousPhase:GaseousPhase");

        surface = surface.withPhases(phase1, phase2);
        CHECK(surface.phases().first.name() == "AqueousPhase");
        CHECK(surface.phases().second.name() == "GaseousPhase");
    }

    WHEN("using constructor Surface(name)")
    {
        Surface surface("AqueousPhase:GaseousPhase");

        CHECK(surface.name() == "AqueousPhase:GaseousPhase");

        surface = surface.withPhases(phase1, phase2);
        CHECK(surface.phases().first.name() == "AqueousPhase");
        CHECK(surface.phases().second.name() == "GaseousPhase");
    }
}
