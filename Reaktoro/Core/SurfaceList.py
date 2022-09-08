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
import pytest

def testSurfaceList():
    #-------------------------------------------------------------------------
    # TESTING CONSTRUCTOR: SurfaceList(formulas)
    #-------------------------------------------------------------------------
    surfaces = SurfaceList([
        Surface("AqueousPhase:GaseousPhase", "AqueousPhase", "GaseousPhase"),
        Surface("Calcite", "Calcite", "Calcite"),
        Surface("Quartz:AqueousPhase", "Quartz", "AqueousPhase"),
    ])

    assert surfaces.size() == 3

    assert surfaces[0].name()  == "AqueousPhase:GaseousPhase"
    assert surfaces[1].name()  == "Calcite"
    assert surfaces[2].name()  == "Quartz:AqueousPhase"

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::find
    #-------------------------------------------------------------------------
    assert surfaces.find("AqueousPhase:GaseousPhase") == 0
    assert surfaces.find("Calcite")                   == 1
    assert surfaces.find("Quartz:AqueousPhase")       == 2

    assert surfaces.find("Xy") >= surfaces.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::findWithName
    #-------------------------------------------------------------------------
    assert surfaces.findWithName("AqueousPhase:GaseousPhase") == 0
    assert surfaces.findWithName("Calcite")                   == 1
    assert surfaces.findWithName("Quartz:AqueousPhase")       == 2

    assert surfaces.findWithName("Xyrium") >= surfaces.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::findWithPhases
    #-------------------------------------------------------------------------
    assert surfaces.findWithPhases("AqueousPhase", "GaseousPhase") == 0
    assert surfaces.findWithPhases("Calcite", "Calcite")           == 1
    assert surfaces.findWithPhases("Quartz", "AqueousPhase")       == 2

    assert surfaces.findWithPhases("Xy", "Zw") >= surfaces.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::index
    #-------------------------------------------------------------------------
    assert surfaces.index("AqueousPhase:GaseousPhase") == 0
    assert surfaces.index("Calcite")                   == 1
    assert surfaces.index("Quartz:AqueousPhase")       == 2

    with pytest.raises(Exception):
        surfaces.index("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::indexWithName
    #-------------------------------------------------------------------------
    assert surfaces.indexWithName("AqueousPhase:GaseousPhase") == 0
    assert surfaces.indexWithName("Calcite")                   == 1
    assert surfaces.indexWithName("Quartz:AqueousPhase")       == 2

    with pytest.raises(Exception):
        surfaces.indexWithName("Xyrium")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::indexWithPhases
    #-------------------------------------------------------------------------
    assert surfaces.indexWithPhases("AqueousPhase", "GaseousPhase") == 0
    assert surfaces.indexWithPhases("Calcite", "Calcite")           == 1
    assert surfaces.indexWithPhases("Quartz", "AqueousPhase")       == 2

    with pytest.raises(Exception):
        surfaces.indexWithPhases("Xy", "Zw")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::get
    #-------------------------------------------------------------------------
    assert surfaces.get("AqueousPhase:GaseousPhase").name() == "AqueousPhase:GaseousPhase"
    assert surfaces.get("Calcite").name()                   == "Calcite"
    assert surfaces.get("Quartz:AqueousPhase").name()       == "Quartz:AqueousPhase"

    with pytest.raises(Exception):
        surfaces.get("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::getWithName
    #-------------------------------------------------------------------------
    assert surfaces.getWithName("AqueousPhase:GaseousPhase").name() == "AqueousPhase:GaseousPhase"
    assert surfaces.getWithName("Calcite").name()                   == "Calcite"
    assert surfaces.getWithName("Quartz:AqueousPhase").name()       == "Quartz:AqueousPhase"

    with pytest.raises(Exception):
        surfaces.getWithName("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::getWithPhases
    #-------------------------------------------------------------------------
    assert surfaces.getWithPhases("AqueousPhase", "GaseousPhase").name() == "AqueousPhase:GaseousPhase"
    assert surfaces.getWithPhases("Calcite", "Calcite").name()           == "Calcite"
    assert surfaces.getWithPhases("Quartz", "AqueousPhase").name()       == "Quartz:AqueousPhase"

    with pytest.raises(Exception):
        surfaces.getWithPhases("Xy", "Zw")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::withNames
    #-------------------------------------------------------------------------
    filtered = surfaces.withNames("Calcite AqueousPhase:GaseousPhase")

    assert filtered.size() == 2
    assert filtered.indexWithName("Calcite") < filtered.size()
    assert filtered.indexWithName("AqueousPhase:GaseousPhase") < filtered.size()

    with pytest.raises(Exception):
        surfaces.withNames("Helium Xylium")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::append
    #-------------------------------------------------------------------------
    surfaces.append(Surface("LiquidPhase:Gel", "LiquidPhase", "Gel"))

    assert surfaces.indexWithName("LiquidPhase:Gel") < surfaces.size()
    assert surfaces.indexWithPhases("LiquidPhase", "Gel") < surfaces.size()
    assert surfaces.indexWithPhases("Gel", "LiquidPhase") < surfaces.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SurfaceList::__iter__
    #-------------------------------------------------------------------------
    for i, surface in enumerate(surfaces):
        assert surface.name() == surfaces[i].name()
