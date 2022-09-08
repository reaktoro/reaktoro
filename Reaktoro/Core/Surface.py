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


def testSurface():
    phase1 = Phase() \
        .withName("AqueousPhase") \
        .withSpecies(SpeciesList("H2O(aq) H+(aq) OH-(aq) CO2(aq) HCO3-(aq) CO3--(aq)")) \
        .withStateOfMatter(StateOfMatter.Liquid)

    phase2 = Phase() \
        .withName("GaseousPhase") \
        .withSpecies(SpeciesList("CO2(g) H2O(g)")) \
        .withStateOfMatter(StateOfMatter.Gas)

    # When using constructor Surface()
    surface = Surface()

    surface = surface.withName("AqueousPhase:GaseousPhase")
    assert surface.name() == "AqueousPhase:GaseousPhase"

    surface = surface.withPhases(phase1, phase2)
    assert surface.phases()[0].name() == "AqueousPhase"
    assert surface.phases()[1].name() == "GaseousPhase"

    # When using constructor Surface(name)
    surface = Surface("AqueousPhase:GaseousPhase")

    assert surface.name() == "AqueousPhase:GaseousPhase"

    surface = surface.withPhases(phase1, phase2)
    assert surface.phases()[0].name() == "AqueousPhase"
    assert surface.phases()[1].name() == "GaseousPhase"
