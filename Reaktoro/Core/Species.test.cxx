// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Singletons/Elements.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Species class", "[Species]")
{
    const auto A = Element().withSymbol("A").withMolarMass(1.0);
    const auto B = Element().withSymbol("B").withMolarMass(2.0);
    const auto C = Element().withSymbol("C").withMolarMass(3.0);
    const auto D = Element().withSymbol("D").withMolarMass(4.0);

    auto species = Species()
        .withName("AB2C3+2(aq)")
        .withFormula("AB2C3+2")
        .withSubstance("AB2C3+2")
        .withElements({{A, 1}, {B, 2}, {C, 3}})
        .withCharge(2.0)
        .withAggregateState(AggregateState::Aqueous)
        .withTags({"tag1", "tag2", "tag3"})
        .withAttachedData(String{"SomeData"})
        ;

    SECTION("Testing attributes of the chemical species")
    {
        REQUIRE( species.name() == "AB2C3+2(aq)" );
        REQUIRE( species.formula() == "AB2C3+2" );
        REQUIRE( species.substance() == "AB2C3+2" );
        REQUIRE( species.elements().size() == 3 );
        REQUIRE( species.elements().coefficient("A") == 1 );
        REQUIRE( species.elements().coefficient("B") == 2 );
        REQUIRE( species.elements().coefficient("C") == 3 );
        REQUIRE( species.charge() == 2.0 );
        REQUIRE( species.tags().size() == 3.0 );
        REQUIRE( species.tags().at(0) == "tag1" );
        REQUIRE( species.tags().at(1) == "tag2" );
        REQUIRE( species.tags().at(2) == "tag3" );
        REQUIRE( species.attachedData().has_value() );
        REQUIRE( species.attachedData().type() == typeid(String) );
        REQUIRE( std::any_cast<String>(species.attachedData()) == "SomeData" );
    }

    SECTION("Testing the standard thermodynamic property functionality of the chemical species")
    {
        const auto T = 300.0;
        const auto P = 1.0e+5;

        REQUIRE_THROWS( species.props(T, P) );

        species = species.withStandardGibbsEnergy(1234.0);

        REQUIRE( species.props(T, P).G0  == 1234.0 );
        REQUIRE( species.props(T, P).H0  == 0.0    );
        REQUIRE( species.props(T, P).V0  == 0.0    );
        REQUIRE( species.props(T, P).Cp0 == 0.0    );
        REQUIRE( species.props(T, P).Cv0 == 0.0    );

        species = species.withStandardThermoModel([](real T, real P) {
            return StandardThermoProps{
                1.0*T*P, // G0
                2.0*T*P, // H0
                3.0*T*P, // V0
                4.0*T*P, // Cp0
                5.0*T*P, // Cv0
            };
        });

        REQUIRE( species.props(T, P).G0  == Approx(1.0*T*P) );
        REQUIRE( species.props(T, P).H0  == Approx(2.0*T*P) );
        REQUIRE( species.props(T, P).V0  == Approx(3.0*T*P) );
        REQUIRE( species.props(T, P).Cp0 == Approx(4.0*T*P) );
        REQUIRE( species.props(T, P).Cv0 == Approx(5.0*T*P) );

        const auto R1 = Species().withName("R1").withStandardGibbsEnergy(0.0);
        const auto R2 = Species().withName("R2").withStandardGibbsEnergy(0.0);

        species = species.withFormationReaction(
            FormationReaction()
                .withReactants({{R1, 1.0}, {R2, 2.0}})
                .withProductStandardVolume(0.1)
                .withReactionThermoModel(
                    [](ReactionThermoProps& props, ReactionThermoArgs args) {
                        const auto& [T, P, dV0] = args;
                        props.dG0 = T + P;
                        props.dH0 = T - P;
                        return props;
                    })
            );

        REQUIRE( species.props(T, P).G0 == species.reaction().createStandardThermoModel()(T, P).G0 );
        REQUIRE( species.props(T, P).H0 == species.reaction().createStandardThermoModel()(T, P).H0 );
    }

    SECTION("Testing automatic construction of chemical species with given chemical formula")
    {
        species = Species("H2O");
        REQUIRE(species.name() == "H2O");
        REQUIRE(species.formula() == "H2O");
        REQUIRE(species.substance() == "H2O");
        REQUIRE(species.charge() == 0);
        REQUIRE(species.molarMass() == Approx(0.01801528));
        REQUIRE(species.aggregateState() == AggregateState::Undefined);
        REQUIRE(species.elements().size() == 2);
        REQUIRE(species.elements().coefficient("H") == 2);
        REQUIRE(species.elements().coefficient("O") == 1);
        REQUIRE(species.tags().empty());

        species = Species("Na+").withName("Na+(aq)").withTags({"aqueous", "cation", "charged"});
        REQUIRE(species.name() == "Na+(aq)");
        REQUIRE(species.formula() == "Na+");
        REQUIRE(species.substance() == "Na+");
        REQUIRE(species.charge() == 1);
        REQUIRE(species.molarMass() == Approx(0.022989769));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 1);
        REQUIRE(species.elements().coefficient("Na") == 1);
        REQUIRE(species.tags().size() == 3);
        REQUIRE(contains(species.tags(), "aqueous"));
        REQUIRE(contains(species.tags(), "cation"));
        REQUIRE(contains(species.tags(), "charged"));

        species = Species("Cl-").withName("Cl-(aq)").withTags({"aqueous", "anion", "charged"});
        REQUIRE(species.name() == "Cl-(aq)");
        REQUIRE(species.formula() == "Cl-");
        REQUIRE(species.substance() == "Cl-");
        REQUIRE(species.charge() == -1);
        REQUIRE(species.molarMass() == Approx(0.035453));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 1);
        REQUIRE(species.elements().coefficient("Cl") == 1);
        REQUIRE(species.tags().size() == 3);
        REQUIRE(contains(species.tags(), "aqueous"));
        REQUIRE(contains(species.tags(), "anion"));
        REQUIRE(contains(species.tags(), "charged"));

        species = Species("CO3--").withName("CO3--(aq)").withTags({"aqueous", "anion", "charged"});
        REQUIRE(species.name() == "CO3--(aq)");
        REQUIRE(species.formula() == "CO3--");
        REQUIRE(species.substance() == "CO3--");
        REQUIRE(species.charge() == -2);
        REQUIRE(species.molarMass() == Approx(0.0600092));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 2);
        REQUIRE(species.elements().coefficient("C") == 1);
        REQUIRE(species.elements().coefficient("O") == 3);
        REQUIRE(species.tags().size() == 3);
        REQUIRE(contains(species.tags(), "aqueous"));
        REQUIRE(contains(species.tags(), "anion"));
        REQUIRE(contains(species.tags(), "charged"));

        species = Species("CaCO3(aq)");
        REQUIRE(species.name() == "CaCO3(aq)");
        REQUIRE(species.formula() == "CaCO3");
        REQUIRE(species.substance() == "CaCO3");
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.charge() == 0);
        REQUIRE(species.molarMass() == Approx(0.1000869));
        REQUIRE(species.elements().size() == 3);
        REQUIRE(species.elements().coefficient("C") == 1);
        REQUIRE(species.elements().coefficient("Ca") == 1);
        REQUIRE(species.elements().coefficient("O") == 3);
        REQUIRE(species.tags().empty());

        species = Species("H+").withName("H+(aq)");
        REQUIRE(species.name() == "H+(aq)");
        REQUIRE(species.formula() == "H+");
        REQUIRE(species.substance() == "H+");
        REQUIRE(species.charge() == 1);
        REQUIRE(species.molarMass() == Approx(0.00100794));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 1);
        REQUIRE(species.elements().coefficient("H") == 1);
        REQUIRE(species.tags().empty());

        species = Species("HCO3-").withTags({"aqueous"});
        REQUIRE(species.name() == "HCO3-");
        REQUIRE(species.formula() == "HCO3-");
        REQUIRE(species.substance() == "HCO3-");
        REQUIRE(species.charge() == -1);
        REQUIRE(species.molarMass() == Approx(0.0610168));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 3);
        REQUIRE(species.elements().coefficient("C") == 1);
        REQUIRE(species.elements().coefficient("H") == 1);
        REQUIRE(species.elements().coefficient("O") == 3);
        REQUIRE(species.tags().size() == 1);
        REQUIRE(contains(species.tags(), "aqueous"));

        species = Species("Fe+++").withTags({"aqueous", "cation", "charged", "iron"});
        REQUIRE(species.name() == "Fe+++");
        REQUIRE(species.formula() == "Fe+++");
        REQUIRE(species.substance() == "Fe+++");
        REQUIRE(species.charge() == 3);
        REQUIRE(species.molarMass() == Approx(0.055847));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 1);
        REQUIRE(species.elements().coefficient("Fe") == 1);
        REQUIRE(species.tags().size() == 4);
        REQUIRE(contains(species.tags(), "aqueous"));
        REQUIRE(contains(species.tags(), "cation"));
        REQUIRE(contains(species.tags(), "charged"));
        REQUIRE(contains(species.tags(), "iron"));
    }

    SECTION("Testing automatic construction of chemical species with given chemical formula containing unknown elements")
    {
        Elements::append(Element().withSymbol("Aa"));
        Elements::append(Element().withSymbol("Bb"));

        species = Species("AaBb2+");
        REQUIRE(species.name() == "AaBb2+");
        REQUIRE(species.formula() == "AaBb2+");
        REQUIRE(species.substance() == "AaBb2+");
        REQUIRE(species.charge() == 1);
        REQUIRE(species.molarMass() == Approx(0.0));
        REQUIRE(species.aggregateState() == AggregateState::Aqueous);
        REQUIRE(species.elements().size() == 2);
        REQUIRE(species.elements().coefficient("Aa") == 1);
        REQUIRE(species.elements().coefficient("Bb") == 2);
        REQUIRE(species.tags().empty());

        REQUIRE_THROWS( Species("RrGgHh") ); // Elements Rr Gg and Hh were not previously appended in the Elements
    }

    SECTION("Testing construction of chemical species with given Species::Attribs")
    {
        const auto T = 300.0;
        const auto P = 1e5;

        WHEN("only formula is provided")
        {
            Species::Attribs attribs;
            attribs.name = "CO3--(aq)";
            attribs.formula = "CO3--";

            species = Species(attribs);

            REQUIRE(species.name() == "CO3--(aq)");
            REQUIRE(species.formula() == "CO3--");
            REQUIRE(species.substance() == "CO3--");
            REQUIRE(species.elements().size() == 2);
            REQUIRE(species.elements().coefficient("C") == 1);
            REQUIRE(species.elements().coefficient("O") == 3);
            REQUIRE(species.charge() == -2.0);
            REQUIRE(species.tags().empty());
            REQUIRE_THROWS(species.props(T, P));
        }

        WHEN("all attributes are given")
        {
            Species::Attribs attribs;
            attribs.name = "CO3--";
            attribs.formula = "CO3--";
            attribs.substance = "CARBONATE";
            attribs.elements = { {Element("C"), 1}, {Element("O"), 3} };
            attribs.charge = -2;
            attribs.aggregate_state = AggregateState::Aqueous;
            attribs.tags = {"aqueous", "solvent"};

            species = Species(attribs);

            REQUIRE(species.name() == "CO3--");
            REQUIRE(species.formula() == "CO3--");
            REQUIRE(species.substance() == "CARBONATE");
            REQUIRE(species.elements().size() == 2);
            REQUIRE(species.elements().coefficient("C") == 1);
            REQUIRE(species.elements().coefficient("O") == 3);
            REQUIRE(species.charge() == -2);
            REQUIRE(species.aggregateState() == AggregateState::Aqueous);
            REQUIRE(species.tags() == attribs.tags);

            WHEN("a StandardThermoModel object is provided for the species")
            {
                attribs.std_thermo_model = [](real T, real P)
                {
                    StandardThermoProps props;
                    props.G0 = 1.234;
                    props.H0 = 2.345;
                    return props;
                };

                species = Species(attribs);

                REQUIRE(species.props(T, P).G0 == 1.234);
                REQUIRE(species.props(T, P).H0 == 2.345);
            }

            WHEN("a FormationReaction object is provided for the species")
            {
                const auto H  = Species("H+").withStandardGibbsEnergy(0.0);
                const auto OH = Species("OH-").withStandardGibbsEnergy(0.0);

                attribs.formation_reaction = FormationReaction()
                    .withReactants({{H, 1.0}, {OH, 1.0}})
                    .withEquilibriumConstant(14.0)
                    .withProductStandardVolume(1e-5)
                    ;

                species = Species(attribs);

                const auto RT = universalGasConstant * T;

                REQUIRE(species.reaction().reactants().size() == 2);
                REQUIRE(species.reaction().reactants().at(0).first.name() == "H+");
                REQUIRE(species.reaction().reactants().at(0).second == 1.0);
                REQUIRE(species.reaction().reactants().at(1).first.name() == "OH-");
                REQUIRE(species.reaction().reactants().at(1).second == 1.0);
                REQUIRE(species.props(T, P).G0 == Approx(-RT*ln10 * 14.0));
                REQUIRE(species.props(T, P).H0 == 0.0);
            }
        }
    }
}
