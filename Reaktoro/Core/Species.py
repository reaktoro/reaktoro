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


def testSpecies():

    # Testing construction of Species with custom elements
    A = Element().withSymbol("A").withMolarMass(1.0)
    B = Element().withSymbol("B").withMolarMass(2.0)
    C = Element().withSymbol("C").withMolarMass(3.0)
    D = Element().withSymbol("D").withMolarMass(4.0)

    species = Species() \
        .withName("AB2C3+2(aq)") \
        .withFormula("AB2C3+2") \
        .withSubstance("AB2C3+2") \
        .withElements(ElementalComposition([(A, 1), (B, 2), (C, 3)])) \
        .withCharge(2.0) \
        .withAggregateState(AggregateState.Aqueous) \
        .withTags("tag1 tag2 tag3")

    assert species.name() == "AB2C3+2(aq)"
    assert species.formula().equivalent("AB2C3+2")
    assert species.substance() == "AB2C3+2"
    assert species.elements().size() == 3
    assert species.elements().coefficient("A") == 1
    assert species.elements().coefficient("B") == 2
    assert species.elements().coefficient("C") == 3
    assert species.charge() == 2.0
    assert species.tags() == ["tag1", "tag2", "tag3"]

    # Testing construction of Species with given chemical formula
    species = Species("CO3--").withName("CO3--(aq)").withTags(["aqueous", "anion", "charged"])
    assert species.name() == "CO3--(aq)"
    assert species.formula().equivalent("CO3--")
    assert species.substance() == "CO3--"
    assert species.charge() == -2
    assert species.molarMass() == pytest.approx(0.0600102972)
    assert species.aggregateState() == AggregateState.Aqueous
    assert species.elements().size() == 2
    assert species.elements().coefficient("C") == 1
    assert species.elements().coefficient("O") == 3
    assert species.tags() == ["aqueous", "anion", "charged"]

    # Testing method Species::props(T, P) with constant standard Gibbs energy
    T = 300.0
    P = 1.0e+5

    with pytest.raises(Exception):
        species.props(T, P)

    species = species.withStandardGibbsEnergy(1234.0)

    assert species.props(T, P).G0[0]  == 1234.0
    assert species.props(T, P).H0[0]  == 0.0
    assert species.props(T, P).V0[0]  == 0.0
    assert species.props(T, P).Cp0[0] == 0.0
    assert species.props(T, P).Cv0[0] == 0.0

    # Testing method Species::props(T, P) with custom standard thermodynamic model
    def model(T, P):
        props = StandardThermoProps()
        props.G0 = 1.0*T*P
        props.H0 = 2.0*T*P
        props.V0 = 3.0*T*P
        props.Cp0 = 4.0*T*P
        props.Cv0 = 5.0*T*P
        return props

    species = species.withStandardThermoModel(StandardThermoModel(model))

    assert species.props(T, P).G0[0]  == pytest.approx(1.0*T*P)
    assert species.props(T, P).H0[0]  == pytest.approx(2.0*T*P)
    assert species.props(T, P).V0[0]  == pytest.approx(3.0*T*P)
    assert species.props(T, P).Cp0[0] == pytest.approx(4.0*T*P)
    assert species.props(T, P).Cv0[0] == pytest.approx(5.0*T*P)

    # Testing method Species::props(T, P) with thermodynamic model of a formation reaction
    R1 = Species().withName("R1").withStandardGibbsEnergy(0.0)
    R2 = Species().withName("R2").withStandardGibbsEnergy(0.0)

    def model(args):
        T, P = args.T, args.P
        props = ReactionThermoProps()
        props.dG0 = T + P
        props.dH0 = T - P
        return props

    species = species.withFormationReaction(
        FormationReaction()
            .withReactants([(R1, 1.0), (R2, 2.0)])
            .withProductStandardVolume(0.1)
            .withReactionThermoModel(ReactionThermoModel(model))
        )

    assert species.props(T, P).G0[0] == species.reaction().createStandardThermoModel()(T, P).G0[0]
    assert species.props(T, P).H0[0] == species.reaction().createStandardThermoModel()(T, P).H0[0]
