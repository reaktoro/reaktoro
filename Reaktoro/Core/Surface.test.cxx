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
    Surface surface;

    surface = surface.withName("AqueousPhase:GaseousPhase");
    CHECK(surface.name() == "AqueousPhase:GaseousPhase");

    surface = surface.withPhaseNames("AqueousPhase", "GaseousPhase");
    CHECK(surface.phaseNames() == Pair<String, String>{"AqueousPhase","GaseousPhase"});

    surface = surface.withPhaseIndices(0, 1);
    CHECK(surface.phaseIndices() == Pair<Index, Index>{0, 1});

    Surface other("SomeName");
    other = other.withPhaseNames("GaseousPhase", "AqueousPhase");
    other = other.withPhaseIndices(1, 0);
    CHECK( surface.equivalent(other) );

    Surface another("SomeName", "GaseousPhase", 1, "AqueousPhase", 0);
    CHECK( surface.equivalent(another) );
}
