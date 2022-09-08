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


def testSurface():

    # When using constructor Surface()
    surface = Surface()

    surface = surface.withName("AqueousPhase:GaseousPhase")
    assert surface.name() == "AqueousPhase:GaseousPhase"

    surface = surface.withPhases("AqueousPhase", "GaseousPhase")
    assert surface.phases()[0] == "AqueousPhase"
    assert surface.phases()[1] == "GaseousPhase"

    other = Surface()
    other = other.withName("SomeName")
    other = other.withPhases("GaseousPhase", "AqueousPhase")

    assert surface.equivalent(other)

    # When using constructor Surface(name)
    surface = Surface("AqueousPhase:GaseousPhase")

    assert surface.name() == "AqueousPhase:GaseousPhase"

    surface = surface.withPhases("AqueousPhase", "GaseousPhase")
    assert surface.phases()[0] == "AqueousPhase"
    assert surface.phases()[1] == "GaseousPhase"
