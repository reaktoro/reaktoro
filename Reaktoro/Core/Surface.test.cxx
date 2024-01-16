// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Surface.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Surface", "[Surface]")
{
    ChemicalProps props;

    auto areafn = [](ChemicalProps const& props) { return 1.23; };

    WHEN("using constructor Surface()")
    {
        Surface surface;

        surface = surface.withName("AqueousPhase:GaseousPhase");
        CHECK( surface.name() == "AqueousPhase:GaseousPhase" );

        surface = surface.withAreaModel(areafn);
        CHECK( surface.areaModel()(props) == 1.23 );
    }

    WHEN("using constructor Surface(name)")
    {
        Surface surface("Quartz");

        CHECK( surface.name() == "Quartz" );
    }

    WHEN("using constructor Surface(name, model)")
    {
        Surface surface("Calcite", areafn);

        CHECK( surface.name() == "Calcite" );
        CHECK( surface.area(props) == 1.23 );
    }
}
