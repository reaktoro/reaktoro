// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Singletons/PeriodicTable.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Species class", "[Species]")
{
    Species species;

    species = Species("H2O");
    REQUIRE(species.formula().equivalent("H2O"));
    REQUIRE(species.name() == "H2O");
    REQUIRE(species.tags().empty());
    REQUIRE(species.molarMass() == Approx(0.01801528));
    REQUIRE(species.charge() == 0);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("H") == 2);
    REQUIRE(species.elementCoefficient("O") == 1);

    species = Species("Na+").withName("Na+(aq)").withTags({"aqueous", "cation", "charged"});
    REQUIRE(species.formula().equivalent("Na+"));
    REQUIRE(species.name() == "Na+(aq)");
    REQUIRE(species.molarMass() == Approx(0.022989769));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Na") == 1);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "cation"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("Cl-").withName("Cl-(aq)").withTags({"aqueous", "anion", "charged"});
    REQUIRE(species.formula().equivalent("Cl-"));
    REQUIRE(species.name() == "Cl-(aq)");
    REQUIRE(species.molarMass() == Approx(0.035453));
    REQUIRE(species.charge() == -1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Cl") == 1);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "anion"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("CO3--").withName("CO3--(aq)").withTags({"aqueous", "anion", "charged"});
    REQUIRE(species.formula().equivalent("CO3-2"));
    REQUIRE(species.name() == "CO3--(aq)");
    REQUIRE(species.molarMass() == Approx(0.0600092));
    REQUIRE(species.charge() == -2);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().size() == 3);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "anion"));
    REQUIRE(contains(species.tags(), "charged"));

    species = Species("CaCO3").withFormula("Ca(CO3)");
    REQUIRE(species.formula().equivalent("Ca(CO3)"));
    REQUIRE(species.name() == "CaCO3");
    REQUIRE(species.molarMass() == Approx(0.1000869));
    REQUIRE(species.charge() == 0);
    REQUIRE(species.elements().size() == 3);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("Ca") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().empty());

    species = Species("H+").withName("H+(aq)");
    REQUIRE(species.formula().equivalent("H+"));
    REQUIRE(species.name() == "H+(aq)");
    REQUIRE(species.molarMass() == Approx(0.00100794));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("H") == 1);
    REQUIRE(species.tags().empty());

    species = Species("HCO3-").withTags({"aqueous"});
    REQUIRE(species.formula().equivalent("HCO3-"));
    REQUIRE(species.name() == "HCO3-");
    REQUIRE(species.molarMass() == Approx(0.0610168));
    REQUIRE(species.charge() == -1);
    REQUIRE(species.elements().size() == 3);
    REQUIRE(species.elementCoefficient("C") == 1);
    REQUIRE(species.elementCoefficient("H") == 1);
    REQUIRE(species.elementCoefficient("O") == 3);
    REQUIRE(species.tags().size() == 1);
    REQUIRE(contains(species.tags(), "aqueous"));

    species = Species("Fe+++").withTags({"aqueous", "cation", "charged", "iron"});
    REQUIRE(species.formula().equivalent("Fe+3"));
    REQUIRE(species.name() == "Fe+++");
    REQUIRE(species.molarMass() == Approx(0.055847));
    REQUIRE(species.charge() == 3);
    REQUIRE(species.elements().size() == 1);
    REQUIRE(species.elementCoefficient("Fe") == 1);
    REQUIRE(species.tags().size() == 4);
    REQUIRE(contains(species.tags(), "aqueous"));
    REQUIRE(contains(species.tags(), "cation"));
    REQUIRE(contains(species.tags(), "charged"));
    REQUIRE(contains(species.tags(), "iron"));

    PeriodicTable::append(Element("Aa"));
    PeriodicTable::append(Element("Bb"));

    species = Species("AaBb2+");
    REQUIRE(species.formula().equivalent("AaBb2+"));
    REQUIRE(species.name() == "AaBb2+");
    REQUIRE(species.molarMass() == Approx(0.0));
    REQUIRE(species.charge() == 1);
    REQUIRE(species.elements().size() == 2);
    REQUIRE(species.elementCoefficient("Aa") == 1);
    REQUIRE(species.elementCoefficient("Bb") == 2);
    REQUIRE(species.tags().empty());

    REQUIRE_THROWS( Species("RrGgHh") ); // Elements Rr Gg and Hh were not previously appended in the PeriodicTable
}
