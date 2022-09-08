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

class ReactionGeneratorUsingClass:
    def __call__(self, phases):
        return [ Reaction().withEquation("H2O(aq) = H+ + OH-") ]

def ReactionGeneratorUsingFunction(phases):
    return [ Reaction().withEquation("H2O(aq) = H2(aq) + 0.5*O2(aq)") ]


def testReactions():
    db = SupcrtDatabase("supcrtbl")

    phases = Phases(db)
    phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- Ca+2 Mg+2 HCO3- CO3-2 CO2(aq) SiO2(aq)") )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhase("Halite") )
    phases.add( MineralPhase("Calcite") )
    phases.add( MineralPhase("Magnesite") )
    phases.add( MineralPhase("Dolomite") )
    phases.add( MineralPhase("Quartz") )

    phases: PhaseList = phases.convert()

    reactions = Reactions()

    reactions.add(db.reaction("Halite = Na+ + Cl-"))
    reactions.add(db.reaction("Calcite"))
    reactions.add(ReactionGeneratorUsingClass())
    reactions.add(ReactionGeneratorUsingFunction)

    converted = reactions.convert(phases)

    assert len(converted) == 4
    assert str(converted[0].equation()) == "Halite = Na+ + Cl-"
    assert str(converted[1].equation()) == "Calcite"
    assert str(converted[2].equation()) == "H2O(aq) = H+ + OH-"
    assert str(converted[3].equation()) == "H2O(aq) = H2(aq) + 0.5*O2(aq)"
