# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *


def generateAreaModel(value: float):
    def ratefn(props: ChemicalProps):
        return value
    return SurfaceAreaModel(ratefn)

class SurfaceGeneratorUsingClass:
    def __call__(self, phases):
        return [ Surface()
            .withName("PhaseA:PhaseB")
            .withAreaModel(generateAreaModel(5.0)) ]

def SurfaceGeneratorUsingFunction(phases):
    return [ Surface()
        .withName("PhaseA:PhaseC")
        .withAreaModel(generateAreaModel(6.0)) ]


def testSurfaces():
    system = ChemicalSystem()
    props = ChemicalProps()

    surface1 = Surface("Surface1")
    surface2 = Surface("Surface2")

    generalsurface1 = GeneralSurface("GeneralSurface1")
    generalsurface2 = GeneralSurface("GeneralSurface2")

    surface1 = surface1.withAreaModel(generateAreaModel(1.0))
    surface2 = surface2.withAreaModel(generateAreaModel(2.0))
    generalsurface1.setAreaModel(generateAreaModel(3.0))
    generalsurface2.setAreaModel(generateAreaModel(4.0))

    surfacesA = Surfaces()
    surfacesA.add(surface1)
    surfacesA.add(surface2)
    surfacesA.add(generalsurface1)
    surfacesA.add(generalsurface2)
    surfacesA.add(SurfaceGeneratorUsingClass())
    surfacesA.add(SurfaceGeneratorUsingFunction)

    surfacesB = Surfaces(
        surface1,
        surface2,
        generalsurface1,
        generalsurface2,
        SurfaceGeneratorUsingClass(),
        SurfaceGeneratorUsingFunction
    )

    def checkSurfacesConversion(surfaces: Surfaces):
        converted = surfaces.convert(system.phases())

        assert len(converted) == 6

        assert converted[0].name() == "Surface1"
        assert converted[1].name() == "Surface2"
        assert converted[2].name() == "GeneralSurface1"
        assert converted[3].name() == "GeneralSurface2"
        assert converted[4].name() == "PhaseA:PhaseB"
        assert converted[5].name() == "PhaseA:PhaseC"

        assert converted[0].area(props).val() == 1.0
        assert converted[1].area(props).val() == 2.0
        assert converted[2].area(props).val() == 3.0
        assert converted[3].area(props).val() == 4.0
        assert converted[4].area(props).val() == 5.0
        assert converted[5].area(props).val() == 6.0

    checkSurfacesConversion(surfacesA)
    checkSurfacesConversion(surfacesB)
