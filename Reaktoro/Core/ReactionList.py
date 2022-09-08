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


def testReactionList():
    species = SpeciesList("H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) NaCl(s) O2(g)")

    db = Database(species)

    reactions = ReactionList()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::append
    #-------------------------------------------------------------------------
    reactions.append(db.reaction("H2O(aq) = H+ + OH-"))
    reactions.append(db.reaction("H2O(aq) = H2(aq) + 0.5*O2(aq)"))
    reactions.append(db.reaction("O2(g) = O2(aq)"))
    reactions.append(db.reaction("NaCl(s) = Na+ + Cl-"))
    reactions.append(db.reaction("NaCl(s)"))

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::size
    #-------------------------------------------------------------------------
    assert reactions.size() == 5

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::operator[](Index)
    #-------------------------------------------------------------------------
    assert reactions[0].name() == "H2O(aq) = H+ + OH-"
    assert reactions[1].name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)"
    assert reactions[2].name() == "O2(g) = O2(aq)"
    assert reactions[3].name() == "NaCl(s) = Na+ + Cl-"
    assert reactions[4].name() == "NaCl(s)"

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::find
    #-------------------------------------------------------------------------
    assert reactions.find("H2O(aq) = H+ + OH-")            == 0
    assert reactions.find("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1
    assert reactions.find("O2(g) = O2(aq)")                == 2
    assert reactions.find("NaCl(s) = Na+ + Cl-")           == 3
    assert reactions.find("NaCl(s)")                       == 4

    assert reactions.find("@#$") >= reactions.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::findWithName
    #-------------------------------------------------------------------------
    assert reactions.findWithName("H2O(aq) = H+ + OH-")            == 0
    assert reactions.findWithName("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1
    assert reactions.findWithName("O2(g) = O2(aq)")                == 2
    assert reactions.findWithName("NaCl(s) = Na+ + Cl-")           == 3
    assert reactions.findWithName("NaCl(s)")                       == 4

    assert reactions.findWithName("@#$") >= reactions.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::index
    #-------------------------------------------------------------------------
    assert reactions.index("H2O(aq) = H+ + OH-")            == 0
    assert reactions.index("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1
    assert reactions.index("O2(g) = O2(aq)")                == 2
    assert reactions.index("NaCl(s) = Na+ + Cl-")           == 3
    assert reactions.index("NaCl(s)")                       == 4

    with pytest.raises(Exception):
        reactions.index("@#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::indexWithName
    #-------------------------------------------------------------------------
    assert reactions.indexWithName("H2O(aq) = H+ + OH-")            == 0
    assert reactions.indexWithName("H2O(aq) = H2(aq) + 0.5*O2(aq)") == 1
    assert reactions.indexWithName("O2(g) = O2(aq)")                == 2
    assert reactions.indexWithName("NaCl(s) = Na+ + Cl-")           == 3
    assert reactions.indexWithName("NaCl(s)")                       == 4

    with pytest.raises(Exception):
        reactions.indexWithName("@#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::withNames
    #-------------------------------------------------------------------------
    filtered = reactions.withNames(["H2O(aq) = H+ + OH-", "H2O(aq) = H2(aq) + 0.5*O2(aq)"])

    assert filtered.size() == 2
    assert filtered[0].name() == "H2O(aq) = H+ + OH-"
    assert filtered[1].name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)"

    with pytest.raises(Exception):
        filtered.withNames("Calcite @#$")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ReactionList::__iter__
    #-------------------------------------------------------------------------
    for i, reaction in enumerate(reactions):
        assert reaction.name() == reactions[i].name()
