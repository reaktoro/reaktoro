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
    REQUIRE(species.symbol() == "H2O");
    REQUIRE(species.name() == "H2O");
    REQUIRE(species.formula() == "H2O");
    REQUIRE(species.tags().empty());
    REQUIRE(species.molarMass() == Approx(0.01801528));
    REQUIRE(species.charge() == 0);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("H") == 2);
    REQUIRE(species.elementCoefficient("O") == 1);

    species = Species("Na+").withSymbol("Na+(aq)").withTags({"aqueous", "cation", "charged"});
    REQUIRE(species.symbol() == "Na+(aq)");
    REQUIRE(species.name() == "Na+");
    REQUIRE(species.formula() == "Na+");
    REQUIRE(species.molarMass() == Approx(0.022989769));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Na") == 1);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "cation"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("Cl-").withSymbol("Cl-(aq)").withTags({"aqueous", "anion", "charged"});
    REQUIRE(species.symbol() == "Cl-(aq)");
    REQUIRE(species.name() == "Cl-");
    REQUIRE(species.formula() == "Cl-");
    REQUIRE(species.molarMass() == Approx(0.035453));
    REQUIRE(species.charge() == -1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Cl") == 1);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "anion"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("CO3--").withSymbol("CO3--(aq)").withTags({"aqueous", "anion", "charged"});
    REQUIRE(species.symbol() == "CO3--(aq)");
    REQUIRE(species.name() == "CO3--");
    REQUIRE(species.formula() == "CO3--");
    REQUIRE(species.molarMass() == Approx(0.0600092));
    REQUIRE(species.charge() == -2);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "anion"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("CaCO3(aq)");
    REQUIRE(species.symbol() == "CaCO3(aq)");
    REQUIRE(species.name() == "CaCO3");
    REQUIRE(species.formula() == "CaCO3");
    REQUIRE(species.molarMass() == Approx(0.1000869));
    REQUIRE(species.charge() == 0);
    REQUIRE(species.elements().size() == 3);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("Ca") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().empty());

    species = Species("H+").withSymbol("H+(aq)");
    REQUIRE(species.symbol() == "H+(aq)");
    REQUIRE(species.name() == "H+");
    REQUIRE(species.formula() == "H+");
    REQUIRE(species.molarMass() == Approx(0.00100794));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("H") == 1);
    REQUIRE(species.tags().empty());

    species = Species("HCO3-").withTags({"aqueous"});
    REQUIRE(species.symbol() == "HCO3-");
    REQUIRE(species.name() == "HCO3-");
    REQUIRE(species.formula() == "HCO3-");
    REQUIRE(species.molarMass() == Approx(0.0610168));
    REQUIRE(species.charge() == -1);
    REQUIRE(species.elements().size() == 3);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("H") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().size() == 1);
    REQUIRE(contains(species.tags(), "aqueous"));

    species = Species("Fe+++").withTags({"aqueous", "cation", "charged", "iron"});
    REQUIRE(species.symbol() == "Fe+++");
    REQUIRE(species.name() == "Fe+++");
    REQUIRE(species.formula() == "Fe+++");
    REQUIRE(species.molarMass() == Approx(0.055847));
    REQUIRE(species.charge() == 3);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Fe") == 1);
    REQUIRE(species.tags().size() == 4);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "cation"));
    REQUIRE(contains(species.tags(), "charged"));
    REQUIRE(contains(species.tags(), "iron"));

    species = Species()
        .withSymbol("H2S(g)")
        .withName("HYDROGEN-SULFIDE")
        .withFormula("H2S")
        .withElementSymbols({{"H", 2}, {"S", 1}})
        .withAggregateState(AggregateState::Gas)
        .withCharge(0.0)
        .withTags({"gaseous"});
    REQUIRE(species.symbol() == "H2S(g)");
    REQUIRE(species.name() == "HYDROGEN-SULFIDE");
    REQUIRE(species.formula() == "H2S");
    REQUIRE(species.molarMass() == Approx(0.0341));
    REQUIRE(species.charge() == 0.0);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("H") == 2);
    REQUIRE(species.elementCoefficient("S") == 1);
    REQUIRE(species.tags().size() == 1);
    REQUIRE(contains(species.tags(), "gaseous"));

    PeriodicTable::append(Element("Aa"));
    PeriodicTable::append(Element("Bb"));

    species = Species("AaBb2+");
    REQUIRE(species.symbol() == "AaBb2+");
    REQUIRE(species.name() == "AaBb2+");
    REQUIRE(species.formula() == "AaBb2+");
    REQUIRE(species.molarMass() == Approx(0.0));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("Aa") == 1);
    REQUIRE(species.elementCoefficient("Bb") == 2);
    REQUIRE(species.tags().empty());

    REQUIRE_THROWS( Species("RrGgHh") ); // Elements Rr Gg and Hh were not previously appended in the PeriodicTable
}
