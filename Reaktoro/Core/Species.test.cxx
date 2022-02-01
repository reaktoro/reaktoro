// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
        CHECK( species.name() == "AB2C3+2(aq)" );
        CHECK( species.formula() == "AB2C3+2" );
        CHECK( species.repr() == "AB2C3+2(aq)" );
        CHECK( species.substance() == "AB2C3+2" );
        CHECK( species.elements().size() == 3 );
        CHECK( species.elements().coefficient("A") == 1 );
        CHECK( species.elements().coefficient("B") == 2 );
        CHECK( species.elements().coefficient("C") == 3 );
        CHECK( species.charge() == 2.0 );
        CHECK( species.tags().size() == 3.0 );
        CHECK( species.tags().at(0) == "tag1" );
        CHECK( species.tags().at(1) == "tag2" );
        CHECK( species.tags().at(2) == "tag3" );
        CHECK( species.attachedData().has_value() );
        CHECK( species.attachedData().type() == typeid(String) );
        CHECK( std::any_cast<String>(species.attachedData()) == "SomeData" );
    }

    SECTION("Testing the standard thermodynamic property functionality of the chemical species")
    {
        const auto T = 300.0;
        const auto P = 1.0e+5;

        CHECK_THROWS( species.standardThermoProps(T, P) );

        species = species.withStandardGibbsEnergy(1234.0);

        CHECK( species.standardThermoProps(T, P).G0  == 1234.0 );
        CHECK( species.standardThermoProps(T, P).H0  == 0.0    );
        CHECK( species.standardThermoProps(T, P).V0  == 0.0    );
        CHECK( species.standardThermoProps(T, P).Cp0 == 0.0    );

        species = species.withStandardThermoModel([](real T, real P) {
            return StandardThermoProps{
                1.0*T*P, // G0
                2.0*T*P, // H0
                3.0*T*P, // V0
                4.0*T*P, // Cp0
                5.0*T*P, // Cv0
            };
        });

        CHECK( species.standardThermoProps(T, P).G0  == Approx(1.0*T*P) );
        CHECK( species.standardThermoProps(T, P).H0  == Approx(2.0*T*P) );
        CHECK( species.standardThermoProps(T, P).V0  == Approx(3.0*T*P) );
        CHECK( species.standardThermoProps(T, P).Cp0 == Approx(4.0*T*P) );

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

        CHECK( species.standardThermoProps(T, P).G0 == species.reaction().createStandardThermoModel()(T, P).G0 );
        CHECK( species.standardThermoProps(T, P).H0 == species.reaction().createStandardThermoModel()(T, P).H0 );
    }

    SECTION("Testing automatic construction of chemical species with given chemical formula")
    {
        species = Species("H2O");
        CHECK(species.name() == "H2O");
        CHECK(species.formula() == "H2O");
        CHECK(species.repr() == "H2O");
        CHECK(species.substance() == "H2O");
        CHECK(species.charge() == 0);
        CHECK(species.molarMass() == Approx(0.01801528));
        CHECK(species.aggregateState() == AggregateState::Undefined);
        CHECK(species.elements().size() == 2);
        CHECK(species.elements().coefficient("H") == 2);
        CHECK(species.elements().coefficient("O") == 1);
        CHECK(species.tags().empty());

        species = Species("CaCO3(s)").withName("Calcite");
        CHECK(species.name() == "Calcite");
        CHECK(species.formula() == "CaCO3");
        CHECK(species.repr() == "Calcite :: CaCO3");
        CHECK(species.substance() == "CaCO3");
        CHECK(species.aggregateState() == AggregateState::Solid);
        CHECK(species.charge() == 0);
        CHECK(species.molarMass() == Approx(0.1000869));
        CHECK(species.elements().size() == 3);
        CHECK(species.elements().coefficient("C") == 1);
        CHECK(species.elements().coefficient("Ca") == 1);
        CHECK(species.elements().coefficient("O") == 3);
        CHECK(species.tags().empty());

        species = Species("Na+").withName("Na+(aq)").withTags({"aqueous", "cation", "charged"});
        CHECK(species.name() == "Na+(aq)");
        CHECK(species.formula() == "Na+");
        CHECK(species.repr() == "Na+(aq)");
        CHECK(species.substance() == "Na+");
        CHECK(species.charge() == 1);
        CHECK(species.molarMass() == Approx(0.0229892194));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 1);
        CHECK(species.elements().coefficient("Na") == 1);
        CHECK(species.tags().size() == 3);
        CHECK(contains(species.tags(), "aqueous"));
        CHECK(contains(species.tags(), "cation"));
        CHECK(contains(species.tags(), "charged"));

        species = Species("Cl-").withName("Cl-(aq)").withTags({"aqueous", "anion", "charged"});
        CHECK(species.name() == "Cl-(aq)");
        CHECK(species.formula() == "Cl-");
        CHECK(species.repr() == "Cl-(aq)");
        CHECK(species.substance() == "Cl-");
        CHECK(species.charge() == -1);
        CHECK(species.molarMass() == Approx(0.035453));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 1);
        CHECK(species.elements().coefficient("Cl") == 1);
        CHECK(species.tags().size() == 3);
        CHECK(contains(species.tags(), "aqueous"));
        CHECK(contains(species.tags(), "anion"));
        CHECK(contains(species.tags(), "charged"));

        species = Species("CO3--").withName("CO3--(aq)").withTags({"aqueous", "anion", "charged"});
        CHECK(species.name() == "CO3--(aq)");
        CHECK(species.formula() == "CO3--");
        CHECK(species.repr() == "CO3--(aq)");
        CHECK(species.substance() == "CO3--");
        CHECK(species.charge() == -2);
        CHECK(species.molarMass() == Approx(0.0600102972));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 2);
        CHECK(species.elements().coefficient("C") == 1);
        CHECK(species.elements().coefficient("O") == 3);
        CHECK(species.tags().size() == 3);
        CHECK(contains(species.tags(), "aqueous"));
        CHECK(contains(species.tags(), "anion"));
        CHECK(contains(species.tags(), "charged"));

        species = Species("CaCO3(aq)");
        CHECK(species.name() == "CaCO3(aq)");
        CHECK(species.formula() == "CaCO3");
        CHECK(species.repr() == "CaCO3(aq)");
        CHECK(species.substance() == "CaCO3");
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.charge() == 0);
        CHECK(species.molarMass() == Approx(0.1000869));
        CHECK(species.elements().size() == 3);
        CHECK(species.elements().coefficient("C") == 1);
        CHECK(species.elements().coefficient("Ca") == 1);
        CHECK(species.elements().coefficient("O") == 3);
        CHECK(species.tags().empty());

        species = Species("H+").withName("H+(aq)");
        CHECK(species.name() == "H+(aq)");
        CHECK(species.formula() == "H+");
        CHECK(species.repr() == "H+(aq)");
        CHECK(species.substance() == "H+");
        CHECK(species.charge() == 1);
        CHECK(species.molarMass() == Approx(0.0010073914));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 1);
        CHECK(species.elements().coefficient("H") == 1);
        CHECK(species.tags().empty());

        species = Species("HCO3-").withTags({"aqueous"});
        CHECK(species.name() == "HCO3-");
        CHECK(species.formula() == "HCO3-");
        CHECK(species.repr() == "HCO3-");
        CHECK(species.substance() == "HCO3-");
        CHECK(species.charge() == -1);
        CHECK(species.molarMass() == Approx(0.0610176886));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 3);
        CHECK(species.elements().coefficient("C") == 1);
        CHECK(species.elements().coefficient("H") == 1);
        CHECK(species.elements().coefficient("O") == 3);
        CHECK(species.tags().size() == 1);
        CHECK(contains(species.tags(), "aqueous"));

        species = Species("Fe+++").withTags({"aqueous", "cation", "charged", "iron"});
        CHECK(species.name() == "Fe+++");
        CHECK(species.formula() == "Fe+++");
        CHECK(species.repr() == "Fe+++");
        CHECK(species.substance() == "Fe+++");
        CHECK(species.charge() == 3);
        CHECK(species.molarMass() == Approx(0.0558453543));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 1);
        CHECK(species.elements().coefficient("Fe") == 1);
        CHECK(species.tags().size() == 4);
        CHECK(contains(species.tags(), "aqueous"));
        CHECK(contains(species.tags(), "cation"));
        CHECK(contains(species.tags(), "charged"));
        CHECK(contains(species.tags(), "iron"));
    }

    SECTION("Testing automatic construction of chemical species with given chemical formula containing unknown elements")
    {
        Elements::append(Element().withSymbol("Aa").withMolarMass(0.01));
        Elements::append(Element().withSymbol("Bb").withMolarMass(0.02));

        species = Species("AaBb2+");
        CHECK(species.name() == "AaBb2+");
        CHECK(species.formula() == "AaBb2+");
        CHECK(species.repr() == "AaBb2+");
        CHECK(species.substance() == "AaBb2+");
        CHECK(species.charge() == 1);
        CHECK(species.molarMass() == Approx(0.0499994514));
        CHECK(species.aggregateState() == AggregateState::Aqueous);
        CHECK(species.elements().size() == 2);
        CHECK(species.elements().coefficient("Aa") == 1);
        CHECK(species.elements().coefficient("Bb") == 2);
        CHECK(species.tags().empty());

        CHECK_THROWS( Species("RrGgHh") ); // Elements Rr Gg and Hh were not previously appended in the Elements
    }

    SECTION("Testing construction of chemical species with given Species::Attribs")
    {
        const auto T = 300.0;
        const auto P = 1e5;

        WHEN("all required attributes are given and a StandardThermoModel object is provided for the species")
        {
            Species::Attribs attribs;
            attribs.name = "CO3--";
            attribs.formula = "CO3--";
            attribs.substance = "CARBONATE";
            attribs.elements = {{Element("C"), 1}, {Element("O"), 3}};
            attribs.charge = -2;
            attribs.aggregate_state = AggregateState::Aqueous;
            attribs.tags = {"aqueous", "carbonate"};
            attribs.std_thermo_model = [](real T, real P)
            {
                StandardThermoProps props;
                props.G0 = 1.234;
                props.H0 = 2.345;
                return props;
            };

            species = Species(attribs);

            CHECK(species.name() == "CO3--");
            CHECK(species.formula() == "CO3--");
            CHECK(species.repr() == "CO3--");
            CHECK(species.substance() == "CARBONATE");
            CHECK(species.elements().size() == 2);
            CHECK(species.elements().coefficient("C") == 1);
            CHECK(species.elements().coefficient("O") == 3);
            CHECK(species.charge() == -2);
            CHECK(species.aggregateState() == AggregateState::Aqueous);
            CHECK(species.tags() == attribs.tags);
            CHECK(species.standardThermoProps(T, P).G0 == 1.234);
            CHECK(species.standardThermoProps(T, P).H0 == 2.345);
        }

        WHEN("all required attributes are given and a FormationReaction object is provided for the species")
        {
            const auto A = Species("H+").withStandardGibbsEnergy(0.0);
            const auto B = Species("OH-").withStandardGibbsEnergy(0.0);

            Species::Attribs attribs;
            attribs.name = "H2O(aq)";
            attribs.formula = "H2O";
            attribs.substance = ""; // empty substance should result in substance equal to formula!
            attribs.elements = {{Element("H"), 2}, {Element("O"), 1}};
            attribs.charge = 0.0;
            attribs.aggregate_state = AggregateState::Aqueous;
            attribs.tags = {"aqueous", "solvent"};
            attribs.formation_reaction = FormationReaction()
                .withReactants({{A, 1.0}, {B, 1.0}})
                .withEquilibriumConstant(14.0)
                .withProductStandardVolume(1e-5);

            species = Species(attribs);

            const auto RT = universalGasConstant * T;

            CHECK(species.name() == "H2O(aq)");
            CHECK(species.formula() == "H2O");
            CHECK(species.repr() == "H2O(aq)");
            CHECK(species.substance() == "H2O");
            CHECK(species.elements().size() == 2);
            CHECK(species.elements().coefficient("H") == 2);
            CHECK(species.elements().coefficient("O") == 1);
            CHECK(species.charge() == 0);
            CHECK(species.aggregateState() == AggregateState::Aqueous);
            CHECK(species.reaction().reactants().size() == 2);
            CHECK(species.reaction().reactants().at(0).first.name() == "H+");
            CHECK(species.reaction().reactants().at(0).second == 1.0);
            CHECK(species.reaction().reactants().at(1).first.name() == "OH-");
            CHECK(species.reaction().reactants().at(1).second == 1.0);
            CHECK(species.standardThermoProps(T, P).G0 == Approx(-RT*ln10 * 14.0));
            CHECK(species.standardThermoProps(T, P).H0 == 0.0);
        }

        WHEN("not all required attributes are given or conflicting attributes are specified")
        {
            using Catch::Contains;

            Species::Attribs attribs;

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::name cannot be empty") );
            attribs.name = "H2O(aq)";

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::formula cannot be empty") );
            attribs.formula = "H2O";

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::elements cannot be empty") );
            attribs.elements = {{Element("H"), 2}, {Element("O"), 1}};

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::aggregate_state cannot be AggregateState::Undefined") );
            attribs.aggregate_state = AggregateState::Aqueous;

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::std_thermo_model and Species::Attribs::formation_reaction cannot be both uninitialized") );
            attribs.std_thermo_model = [](real T, real P)
            {
                StandardThermoProps props;
                props.G0 = 1.234;
                props.H0 = 2.345;
                return props;
            };

            CHECK_NOTHROW( Species(attribs) );

            const auto A = Species("H+").withStandardGibbsEnergy(0.0);
            const auto B = Species("OH-").withStandardGibbsEnergy(0.0);
            attribs.formation_reaction = FormationReaction()
                .withReactants({{A, 1.0}, {B, 1.0}})
                .withEquilibriumConstant(14.0)
                .withProductStandardVolume(1e-5);

            CHECK_THROWS_WITH( Species(attribs), Contains("Species::Attribs::std_thermo_model and Species::Attribs::formation_reaction cannot be both initialized") );
        }
    }
}
