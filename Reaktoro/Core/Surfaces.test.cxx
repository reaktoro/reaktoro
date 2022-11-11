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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Surfaces.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

namespace test {

auto generateAreaModel(double value)
{
    return [=](ChemicalProps const& props) { return value; };
}

class SurfaceGeneratorUsingClass
{
public:
    auto operator()(PhaseList const& phases) const -> Vec<Surface>
    {
        return { Surface()
            .withName("PhaseA:PhaseB")
            .withAreaModel(generateAreaModel(5.0))
        };
    }
};

auto SurfaceGeneratorUsingFunction(PhaseList const& phases) -> Vec<Surface>
{
    return { Surface()
        .withName("PhaseA:PhaseC")
        .withAreaModel(generateAreaModel(6.0))
    };
}

} // namespace test

TEST_CASE("Testing Surfaces class", "[Surfaces]")
{
    ChemicalSystem system;
    ChemicalProps props;

    Surface surface1("Surface1");
    Surface surface2("Surface2");

    GeneralSurface generalsurface1("GeneralSurface1");
    GeneralSurface generalsurface2("GeneralSurface2");

    surface1 = surface1.withAreaModel(test::generateAreaModel(1.0));
    surface2 = surface2.withAreaModel(test::generateAreaModel(2.0));
    generalsurface1.setAreaModel(test::generateAreaModel(3.0));
    generalsurface2.setAreaModel(test::generateAreaModel(4.0));

    Surfaces surfacesA;
    surfacesA.add(surface1);
    surfacesA.add(surface2);
    surfacesA.add(generalsurface1);
    surfacesA.add(generalsurface2);
    surfacesA.add(test::SurfaceGeneratorUsingClass());
    surfacesA.add(test::SurfaceGeneratorUsingFunction);

    Surfaces surfacesB(
        surface1,
        surface2,
        generalsurface1,
        generalsurface2,
        test::SurfaceGeneratorUsingClass(),
        test::SurfaceGeneratorUsingFunction
    );

    auto checkSurfacesConversion = [&](Surfaces const& surfaces)
    {
        auto converted = surfaces.convert(system.phases());

        CHECK( converted.size() == 6 );

        CHECK( converted[0].name() == "Surface1" );
        CHECK( converted[1].name() == "Surface2" );
        CHECK( converted[2].name() == "GeneralSurface1" );
        CHECK( converted[3].name() == "GeneralSurface2" );
        CHECK( converted[4].name() == "PhaseA:PhaseB" );
        CHECK( converted[5].name() == "PhaseA:PhaseC" );

        CHECK( converted[0].area(props).val() == 1.0 );
        CHECK( converted[1].area(props).val() == 2.0 );
        CHECK( converted[2].area(props).val() == 3.0 );
        CHECK( converted[3].area(props).val() == 4.0 );
        CHECK( converted[4].area(props).val() == 5.0 );
        CHECK( converted[5].area(props).val() == 6.0 );
    };

    checkSurfacesConversion(surfacesA);
    checkSurfacesConversion(surfacesB);
}
