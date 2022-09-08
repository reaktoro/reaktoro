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


def testSurfaces():
    surfaces = Surfaces()
    surfaces.add("Calcite", "AqueousPhase")
    surfaces.add("GaseousPhase", "AqueousPhase")
    surfaces.add("Quartz")

    assert len(surfaces.data()) == 3

    assert surfaces.data()[0] == ("Calcite", "AqueousPhase")
    assert surfaces.data()[1] == ("GaseousPhase", "AqueousPhase")
    assert surfaces.data()[2] == ("Quartz", "Quartz")

    phases = PhaseList([
        Phase().withName("AqueousPhase"), # # 0
        Phase().withName("GaseousPhase"), # # 1
        Phase().withName("LiquidPhase"),  # # 2
        Phase().withName("Calcite"),      # # 3
        Phase().withName("Quartz"),       # # 4
        Phase().withName("Dolomite")      # # 5
    ])

    converted = surfaces.convert(phases)

    assert len(converted) == 3

    assert converted[0].name() == "Calcite:AqueousPhase"
    assert converted[1].name() == "GaseousPhase:AqueousPhase"
    assert converted[2].name() == "Quartz"

    assert converted[0].phaseNames() == ("Calcite", "AqueousPhase")
    assert converted[1].phaseNames() == ("GaseousPhase", "AqueousPhase")
    assert converted[2].phaseNames() == ("Quartz", "Quartz")

    assert converted[0].phaseIndices() == (3, 0)
    assert converted[1].phaseIndices() == (1, 0)
    assert converted[2].phaseIndices() == (4, 4)
