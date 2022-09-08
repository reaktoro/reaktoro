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
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/SurfaceList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SurfaceList", "[SurfaceList]")
{
    SurfaceList surfaces;
    SurfaceList filtered;

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: SurfaceList(formulas)
    //-------------------------------------------------------------------------
    surfaces = SurfaceList({
        Surface("AqueousPhase:GaseousPhase").withPhases("AqueousPhase", "GaseousPhase"),
        Surface("Calcite").withPhases("Calcite", "Calcite"),
        Surface("Quartz:AqueousPhase").withPhases("Quartz", "AqueousPhase"),
    });

    REQUIRE( surfaces.size() == 3 );

    REQUIRE( surfaces[0].name() == "AqueousPhase:GaseousPhase" );
    REQUIRE( surfaces[1].name() == "Calcite" );
    REQUIRE( surfaces[2].name() == "Quartz:AqueousPhase" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::find
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.find("AqueousPhase:GaseousPhase") == 0 );
    REQUIRE( surfaces.find("Calcite")                   == 1 );
    REQUIRE( surfaces.find("Quartz:AqueousPhase")       == 2 );

    REQUIRE( surfaces.find("Xy") >= surfaces.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.findWithName("AqueousPhase:GaseousPhase") == 0 );
    REQUIRE( surfaces.findWithName("Calcite")                   == 1 );
    REQUIRE( surfaces.findWithName("Quartz:AqueousPhase")       == 2 );

    REQUIRE( surfaces.findWithName("Xyrium") >= surfaces.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::findWithPhases
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.findWithPhases("AqueousPhase", "GaseousPhase") == 0 );
    REQUIRE( surfaces.findWithPhases("Calcite", "Calcite")           == 1 );
    REQUIRE( surfaces.findWithPhases("Quartz", "AqueousPhase")       == 2 );

    REQUIRE( surfaces.findWithPhases("Xy", "Zw") >= surfaces.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::index
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.index("AqueousPhase:GaseousPhase") == 0 );
    REQUIRE( surfaces.index("Calcite")                   == 1 );
    REQUIRE( surfaces.index("Quartz:AqueousPhase")       == 2 );

    REQUIRE_THROWS( surfaces.index("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::indexWithName
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.indexWithName("AqueousPhase:GaseousPhase") == 0 );
    REQUIRE( surfaces.indexWithName("Calcite")                   == 1 );
    REQUIRE( surfaces.indexWithName("Quartz:AqueousPhase")       == 2 );

    REQUIRE_THROWS( surfaces.indexWithName("Xyrium") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::indexWithPhases
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.indexWithPhases("AqueousPhase", "GaseousPhase") == 0 );
    REQUIRE( surfaces.indexWithPhases("Calcite", "Calcite")           == 1 );
    REQUIRE( surfaces.indexWithPhases("Quartz", "AqueousPhase")       == 2 );

    REQUIRE_THROWS( surfaces.indexWithPhases("Xy", "Zw") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::get
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.get("AqueousPhase:GaseousPhase").name() == "AqueousPhase:GaseousPhase" );
    REQUIRE( surfaces.get("Calcite").name()                   == "Calcite" );
    REQUIRE( surfaces.get("Quartz:AqueousPhase").name()       == "Quartz:AqueousPhase" );

    REQUIRE_THROWS( surfaces.get("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::getWithName
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.getWithName("AqueousPhase:GaseousPhase").name() == "AqueousPhase:GaseousPhase" );
    REQUIRE( surfaces.getWithName("Calcite").name()                   == "Calcite" );
    REQUIRE( surfaces.getWithName("Quartz:AqueousPhase").name()       == "Quartz:AqueousPhase" );

    REQUIRE_THROWS( surfaces.getWithName("Xy") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::getWithPhases
    //-------------------------------------------------------------------------
    REQUIRE( surfaces.getWithPhases("AqueousPhase", "GaseousPhase").name() == "AqueousPhase:GaseousPhase" );
    REQUIRE( surfaces.getWithPhases("Calcite", "Calcite").name()           == "Calcite" );
    REQUIRE( surfaces.getWithPhases("Quartz", "AqueousPhase").name()       == "Quartz:AqueousPhase" );

    REQUIRE_THROWS( surfaces.getWithPhases("Xy", "Zw") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::withNames
    //-------------------------------------------------------------------------
    filtered = surfaces.withNames("Calcite AqueousPhase:GaseousPhase");

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered.indexWithName("Calcite") < filtered.size() );
    REQUIRE( filtered.indexWithName("AqueousPhase:GaseousPhase") < filtered.size() );

    REQUIRE_THROWS( surfaces.withNames("Helium Xylium") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::append
    //-------------------------------------------------------------------------
    surfaces.append(Surface("LiquidPhase:Gel").withPhases("LiquidPhase", "Gel"));

    REQUIRE( surfaces.indexWithName("LiquidPhase:Gel") < surfaces.size() );
    REQUIRE( surfaces.indexWithPhases("LiquidPhase", "Gel") < surfaces.size() );
    REQUIRE( surfaces.indexWithPhases("Gel", "LiquidPhase") < surfaces.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SurfaceList::begin|end
    //-------------------------------------------------------------------------
    for(auto [i, surface] : enumerate(surfaces))
        REQUIRE( surface.name() == surfaces[i].name() );
}
