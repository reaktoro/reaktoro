# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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


def createReactionRateModelAsFunction(value: float):
    def ratefn(props: ChemicalProps):
        print("ratefn1")
        return ReactionRate(value)
    return ratefn


def createReactionRateModel(value: float) -> ReactionRateModel:
    return ReactionRateModel(createReactionRateModelAsFunction(value))


def createReactionRateModelGeneratorAsFunction(value: float):
    def genfn(args: ReactionRateModelGeneratorArgs):
        def ratefn(props: ChemicalProps):
            print("ratefn3")
            return value
        return ratefn
    return genfn


class ReactionGeneratorUsingClass:
    def __call__(self, phases):
        return [ Reaction()
            .withEquation("H2O(aq) = H+ + OH-")
            .withRateModel(createReactionRateModel(6.0)) ]


def ReactionGeneratorUsingFunction(phases):
    return [ Reaction()
        .withEquation("H2O(aq) = H2(aq) + 0.5*O2(aq)")
        .withRateModel(createReactionRateModel(7.0)) ]


def testGeneralReaction():

    species = SpeciesList("CO2(aq) CO2(g)")

    db = Database(species)

    props = ChemicalProps()

    args = ReactionGeneratorArgs(db, species, PhaseList(), SurfaceList())

    reaction = GeneralReaction("CO2(g) = CO2(aq)")

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Check when a Python callable object with single argument of type ChemicalProps is provided to setRateModel (which will be converted to a ReactionRateModel object)
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    try: reaction.setRateModel(createReactionRateModelAsFunction(1.0))
    except: pytest.fail("Expecting GeneralReaction.setRateModel to work with a callable with single argument of type ChemicalProps")

    rxn = reaction.convert(args)
    assert rxn.rate(props).val() == 1.0

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Check when a ReactionRateModel object is provided to setRateModel
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    try: reaction.setRateModel(createReactionRateModel(2.0))
    except: pytest.fail("Expecting GeneralReaction.setRateModel to work with a ReactionRateModel object")

    rxn = reaction.convert(args)
    assert rxn.rate(props).val() == 2.0

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Check when a Python callable object with single argument of type ReactionRateModelGeneratorArgs is provided to setRateModel (which will be converted to a ReactionRateModelGenerator object)
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    try: reaction.setRateModel(createReactionRateModelGeneratorAsFunction(3.0))
    except: pytest.fail("Expecting GeneralReaction.setRateModel to work with a callable with single argument of type ReactionRateModelGeneratorArgs")

    rxn = reaction.convert(args)
    assert rxn.rate(props).val() == 3.0


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

    phases = PhaseList(phases.convert())

    system = ChemicalSystem(db, phases)
    props = ChemicalProps(system)

    reaction1 = db.reaction("Halite = Na+ + Cl-")
    reaction2 = db.reaction("Calcite")

    reaction1 = reaction1.withRateModel(createReactionRateModel(1.0))
    reaction2 = reaction2.withRateModel(createReactionRateModel(2.0))

    generalreaction1 = GeneralReaction("CO2(g) = CO2(aq)").setRateModel(createReactionRateModel(3.0))
    generalreaction2 = GeneralReaction("HCO3- + H+ = CO2(aq) + H2O(aq)").setRateModel(createReactionRateModelAsFunction(4.0))
    generalreaction3 = GeneralReaction("SiO2(aq) = Quartz").setRateModel(createReactionRateModelGeneratorAsFunction(5.0))

    reactionsA = Reactions()
    reactionsA.add(reaction1)
    reactionsA.add(reaction2)
    reactionsA.add(generalreaction1)
    reactionsA.add(generalreaction2)
    reactionsA.add(generalreaction3)
    reactionsA.add(ReactionGeneratorUsingClass())
    reactionsA.add(ReactionGeneratorUsingFunction)

    reactionsB = Reactions(
        reaction1,
        reaction2,
        generalreaction1,
        generalreaction2,
        generalreaction3,
        ReactionGeneratorUsingClass(),
        ReactionGeneratorUsingFunction
    )

    def checkReactionsConversion(reactions: Reactions):
        args = ReactionGeneratorArgs(
            system.database(),
            system.species(),
            system.phases(),
            system.surfaces()
        )

        converted: list[Reaction] = reactions.convert(args)

        assert len(converted) == 7

        assert str(converted[0].equation()) == "Halite = Na+ + Cl-"
        assert str(converted[1].equation()) == "Calcite"
        assert str(converted[2].equation()) == "CO2(g) = CO2(aq)"
        assert str(converted[3].equation()) == "HCO3- + H+ = CO2(aq) + H2O(aq)"
        assert str(converted[4].equation()) == "SiO2(aq) = Quartz"
        assert str(converted[5].equation()) == "H2O(aq) = H+ + OH-"
        assert str(converted[6].equation()) == "H2O(aq) = H2(aq) + 0.5*O2(aq)"

        assert converted[0].rate(props).val() == 1.0
        assert converted[1].rate(props).val() == 2.0
        assert converted[2].rate(props).val() == 3.0
        assert converted[3].rate(props).val() == 4.0
        assert converted[4].rate(props).val() == 5.0
        assert converted[5].rate(props).val() == 6.0
        assert converted[6].rate(props).val() == 7.0

    checkReactionsConversion(reactionsA)
    checkReactionsConversion(reactionsB)
