# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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


def testPhaseList():
    phases = PhaseList()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::append
    #-------------------------------------------------------------------------
    phases.append(Phase()
        .withName("AqueousPhase")
        .withSpecies(SpeciesList("H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq)"))
        .withStateOfMatter(StateOfMatter.Liquid))

    phases.append(Phase()
        .withName("GaseousPhase")
        .withSpecies(SpeciesList("H2O(g) CO2(g) CH4(g) O2(g) H2(g) CO(g)"))
        .withStateOfMatter(StateOfMatter.Gas))

    phases.append(Phase()
        .withName("Calcite")
        .withSpecies([ Species("CaCO3(s)").withName("Calcite") ])
        .withStateOfMatter(StateOfMatter.Solid))

    phases.append(Phase()
        .withName("Halite")
        .withSpecies([ Species("NaCl(s)").withName("Halite") ])
        .withStateOfMatter(StateOfMatter.Solid))

    phases.append(Phase()
        .withName("Magnesite")
        .withSpecies([ Species("MgCO3(s)").withName("Magnesite") ])
        .withStateOfMatter(StateOfMatter.Solid))

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::size
    #-------------------------------------------------------------------------
    assert phases.size() == 5

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::operator[](Index)
    #-------------------------------------------------------------------------
    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Calcite"
    assert phases[3].name() == "Halite"
    assert phases[4].name() == "Magnesite"

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::find
    #-------------------------------------------------------------------------
    assert phases.find("AqueousPhase") == 0
    assert phases.find("GaseousPhase") == 1
    assert phases.find("Calcite")      == 2
    assert phases.find("Halite")       == 3
    assert phases.find("Magnesite")    == 4

    assert phases.find("@#$") >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::findWithName
    #-------------------------------------------------------------------------
    assert phases.findWithName("AqueousPhase") == 0
    assert phases.findWithName("GaseousPhase") == 1
    assert phases.findWithName("Calcite")      == 2
    assert phases.findWithName("Halite")       == 3
    assert phases.findWithName("Magnesite")    == 4

    assert phases.findWithName("@#$") >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::findWithSpecies(index)
    #-------------------------------------------------------------------------
    assert phases.findWithSpecies(0)  == 0  #  idx(0) -- H2O(aq)
    assert phases.findWithSpecies(5)  == 0  #  idx(5) -- Na+
    assert phases.findWithSpecies(8)  == 1  #  idx(8) -- H2O(g)
    assert phases.findWithSpecies(13) == 1  # idx(13) -- CO(g)
    assert phases.findWithSpecies(14) == 2  # idx(14) -- Calcite
    assert phases.findWithSpecies(15) == 3  # idx(15) -- Halite
    assert phases.findWithSpecies(16) == 4  # idx(16) -- Magnesite

    assert phases.findWithSpecies(100) >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::findWithSpecies(name)
    #-------------------------------------------------------------------------
    assert phases.findWithSpecies("H2O(aq)")   == 0
    assert phases.findWithSpecies("Na+")       == 0
    assert phases.findWithSpecies("H2O(g)")    == 1
    assert phases.findWithSpecies("CO(g)")     == 1
    assert phases.findWithSpecies("Calcite")   == 2
    assert phases.findWithSpecies("Halite")    == 3
    assert phases.findWithSpecies("Magnesite") == 4

    assert phases.findWithSpecies("@#$") >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::findWithAggregateState
    #-------------------------------------------------------------------------
    assert phases.findWithAggregateState(AggregateState.Aqueous) == 0
    assert phases.findWithAggregateState(AggregateState.Gas)     == 1
    assert phases.findWithAggregateState(AggregateState.Solid)   == 2

    assert phases.findWithAggregateState(AggregateState.Undefined) >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::findWithStateOfMatter
    #-------------------------------------------------------------------------
    assert phases.findWithStateOfMatter(StateOfMatter.Liquid) == 0
    assert phases.findWithStateOfMatter(StateOfMatter.Gas)    == 1
    assert phases.findWithStateOfMatter(StateOfMatter.Solid)  == 2

    assert phases.findWithStateOfMatter(StateOfMatter.Plasma) >= phases.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::index
    #-------------------------------------------------------------------------
    assert phases.index("AqueousPhase") == 0
    assert phases.index("GaseousPhase") == 1
    assert phases.index("Calcite")      == 2
    assert phases.index("Halite")       == 3
    assert phases.index("Magnesite")    == 4

    with pytest.raises(Exception):
        phases.index("@#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::indexWithName
    #-------------------------------------------------------------------------
    assert phases.indexWithName("AqueousPhase") == 0
    assert phases.indexWithName("GaseousPhase") == 1
    assert phases.indexWithName("Calcite")      == 2
    assert phases.indexWithName("Halite")       == 3
    assert phases.indexWithName("Magnesite")    == 4

    with pytest.raises(Exception):
        phases.indexWithName("@#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::indexWithSpecies(index)
    #-------------------------------------------------------------------------
    assert phases.indexWithSpecies(0)  == 0  #  idx(0) -- H2O(aq)
    assert phases.indexWithSpecies(5)  == 0  #  idx(5) -- Na+
    assert phases.indexWithSpecies(8)  == 1  #  idx(8) -- H2O(g)
    assert phases.indexWithSpecies(13) == 1  # idx(13) -- CO(g)
    assert phases.indexWithSpecies(14) == 2  # idx(14) -- Calcite
    assert phases.indexWithSpecies(15) == 3  # idx(15) -- Halite
    assert phases.indexWithSpecies(16) == 4  # idx(16) -- Magnesite

    with pytest.raises(Exception):
        phases.indexWithSpecies(100)

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::indexWithSpecies(name)
    #-------------------------------------------------------------------------
    assert phases.indexWithSpecies("H2O(aq)")   == 0
    assert phases.indexWithSpecies("Na+")       == 0
    assert phases.indexWithSpecies("H2O(g)")    == 1
    assert phases.indexWithSpecies("CO(g)")     == 1
    assert phases.indexWithSpecies("Calcite")   == 2
    assert phases.indexWithSpecies("Halite")    == 3
    assert phases.indexWithSpecies("Magnesite") == 4

    with pytest.raises(Exception):
        phases.indexWithSpecies("@#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::indexWithAggregateState
    #-------------------------------------------------------------------------
    assert phases.indexWithAggregateState(AggregateState.Aqueous) == 0
    assert phases.indexWithAggregateState(AggregateState.Gas)     == 1
    assert phases.indexWithAggregateState(AggregateState.Solid)   == 2

    with pytest.raises(Exception):
        phases.indexWithAggregateState(AggregateState.Undefined)

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::indexWithStateOfMatter
    #-------------------------------------------------------------------------
    assert phases.indexWithStateOfMatter(StateOfMatter.Liquid) == 0
    assert phases.indexWithStateOfMatter(StateOfMatter.Gas)    == 1
    assert phases.indexWithStateOfMatter(StateOfMatter.Solid)  == 2

    with pytest.raises(Exception):
        phases.indexWithStateOfMatter(StateOfMatter.Plasma)

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::withNames
    #-------------------------------------------------------------------------
    filtered = phases.withNames("AqueousPhase GaseousPhase")

    assert filtered.size() == 2
    assert filtered[0].name() == "AqueousPhase"
    assert filtered[1].name() == "GaseousPhase"

    with pytest.raises(Exception):
        filtered.withNames("Calcite @#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::withStateOfMatter
    #-------------------------------------------------------------------------
    filtered = phases.withStateOfMatter(StateOfMatter.Solid)

    assert filtered.size() == 3
    assert filtered[0].name() == "Calcite"
    assert filtered[1].name() == "Halite"
    assert filtered[2].name() == "Magnesite"

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::withAggregateState
    #-------------------------------------------------------------------------
    filtered = phases.withAggregateState(AggregateState.Aqueous)

    assert filtered.size() == 1
    assert filtered[0].name() == "AqueousPhase"

    filtered = phases.withAggregateState(AggregateState.Solid)

    assert filtered.size() == 3
    assert filtered[0].name() == "Calcite"
    assert filtered[1].name() == "Halite"
    assert filtered[2].name() == "Magnesite"

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::numSpeciesUntilPhase
    #-------------------------------------------------------------------------
    assert phases.numSpeciesUntilPhase(0) == 0
    assert phases.numSpeciesUntilPhase(1) == 8
    assert phases.numSpeciesUntilPhase(2) == 8 + 6
    assert phases.numSpeciesUntilPhase(3) == 8 + 6 + 1
    assert phases.numSpeciesUntilPhase(4) == 8 + 6 + 1 + 1
    assert phases.numSpeciesUntilPhase(5) == 8 + 6 + 1 + 1 + 1

    #-------------------------------------------------------------------------
    # TESTING METHOD: PhaseList::__iter__
    #-------------------------------------------------------------------------
    for i, phase in enumerate(phases):
        assert phase.name() == phases[i].name()
