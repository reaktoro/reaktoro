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
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Surface.hpp>
#include <Reaktoro/Core/Surfaces.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Surfaces", "[Surfaces]")
{
    Surfaces surfaces;
    surfaces.add("Calcite", "AqueousPhase");
    surfaces.add("GaseousPhase", "AqueousPhase");
    surfaces.add("Quartz");

    CHECK( surfaces.data().size() == 3 );

    CHECK( surfaces.data()[0] == Pair<String, String>{"Calcite", "AqueousPhase"} );
    CHECK( surfaces.data()[1] == Pair<String, String>{"GaseousPhase", "AqueousPhase"} );
    CHECK( surfaces.data()[2] == Pair<String, String>{"Quartz", "Quartz"} );

    auto const phases = PhaseList({
        Phase().withName("AqueousPhase"), // # 0
        Phase().withName("GaseousPhase"), // # 1
        Phase().withName("LiquidPhase"),  // # 2
        Phase().withName("Calcite"),      // # 3
        Phase().withName("Quartz"),       // # 4
        Phase().withName("Dolomite")      // # 5
    });

    auto const converted = surfaces.convert(phases);

    CHECK( converted.size() == 3 );

    CHECK( converted[0].name() == "Calcite:AqueousPhase" );
    CHECK( converted[1].name() == "GaseousPhase:AqueousPhase" );
    CHECK( converted[2].name() == "Quartz" );

    CHECK( converted[0].phaseNames() == Pair<String, String>{"Calcite", "AqueousPhase"} );
    CHECK( converted[1].phaseNames() == Pair<String, String>{"GaseousPhase", "AqueousPhase"} );
    CHECK( converted[2].phaseNames() == Pair<String, String>{"Quartz", "Quartz"} );

    CHECK( converted[0].phaseIndices() == Pair<Index, Index>{3, 0} );
    CHECK( converted[1].phaseIndices() == Pair<Index, Index>{1, 0} );
    CHECK( converted[2].phaseIndices() == Pair<Index, Index>{4, 4} );
}
