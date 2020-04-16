// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Singletons/PeriodicTable.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Species class", "[Species]")
{
    Species species;

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

    species = Species()
        .withName("H2S(g)")
        .withSubstance("HYDROGEN-SULFIDE")
        .withFormula("H2S")
        .withElements(Map<String, double>{ {"H", 2}, {"S", 1} })
        .withAggregateState(AggregateState::Gas)
        .withCharge(0.0)
        .withTags({"gaseous"});
    REQUIRE(species.name() == "H2S(g)");
    REQUIRE(species.formula() == "H2S");
    REQUIRE(species.substance() == "HYDROGEN-SULFIDE");
    REQUIRE(species.charge() == 0.0);
    REQUIRE(species.molarMass() == Approx(0.03408188));
    REQUIRE(species.aggregateState() == AggregateState::Gas);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elements().coefficient("H") == 2);
    REQUIRE(species.elements().coefficient("S") == 1);
    REQUIRE(species.tags().size() == 1);
    REQUIRE(contains(species.tags(), "gaseous"));

    PeriodicTable::append(Element("Aa"));
    PeriodicTable::append(Element("Bb"));

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

    REQUIRE_THROWS( Species("RrGgHh") ); // Elements Rr Gg and Hh were not previously appended in the PeriodicTable
}
