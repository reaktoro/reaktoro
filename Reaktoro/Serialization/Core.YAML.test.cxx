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
using namespace Catch;

// Reaktoro includes
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Serialization/Common.YAML.hpp>
#include <Reaktoro/Serialization/Core.YAML.hpp>

using namespace Reaktoro;

TEST_CASE("Testing YAML encoder/decoder for AggregateState", "[Core.yaml]")
{
    yaml node;
    AggregateState aggstate;

    node = AggregateState::Aqueous;
    CHECK( node.repr() == "Aqueous" );
    aggstate = node;
    CHECK( aggstate == AggregateState::Aqueous );

    node = AggregateState::Gas;
    CHECK( node.repr() == "Gas" );
    aggstate = node;
    CHECK( aggstate == AggregateState::Gas );

    node = AggregateState::Solid;
    CHECK( node.repr() == "Solid" );
    aggstate = node;
    CHECK( aggstate == AggregateState::Solid );

    node = AggregateState::Undefined;
    CHECK( node.repr() == "Undefined" );
    aggstate = node;
    CHECK( aggstate == AggregateState::Undefined );
}

TEST_CASE("Testing YAML encoder/decoder for ChemicalFormula", "[Core.yaml]")
{
    yaml node;
    ChemicalFormula formula;

    node = ChemicalFormula("H2O");
    CHECK( node.repr() == "H2O" );
    formula = node;
    CHECK( formula.equivalent("H2O") );

    node = ChemicalFormula("Ca++");
    CHECK( node.repr() == "Ca++" );
    formula = node;
    CHECK( formula.equivalent("Ca++") );
}

TEST_CASE("Testing YAML encoder/decoder for ChemicalSystem", "[Core.yaml]")
{
    yaml node;
    ChemicalSystem chemicalsystem;

    // TODO: Implement tests for YAML encoding/decoding of ChemicalSystem objects.
}

TEST_CASE("Testing YAML encoder/decoder for Database", "[Core.yaml]")
{
    auto contents = R"(
Elements:
  - Symbol: A
    MolarMass: 1.0
  - Symbol: B
    MolarMass: 2.0
  - Symbol: C
    MolarMass: 3.0
Species:
  - Name: A2B(aq)
    Formula: A2B
    Elements: 2:A 1:B
    AggregateState: Aqueous
    FormationReaction:
      Reactants: 2:A(aq) 1:B(aq)
      ReactionThermoModel:
        ConstLgK:
          lgKr: 5.0
  - Name: A(aq)
    Formula: A
    Elements: 1:A
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 1.0 }
  - Name: B(aq)
    Formula: B
    Elements: 1:B
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 2.0 }
  - Name: C(aq)
    Formula: C
    Elements: 1:C
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 3.0 }
  - Name: A3B5C3(aq)
    Formula: A3B5C3
    Elements: 3:A 5:B 3:C
    AggregateState: Aqueous
    FormationReaction:
      Reactants: 1:A2B(aq) 1:A(aq) 4:B(aq) 3:C(aq)
      ReactionThermoModel:
        ConstLgK:
          lgKr: 7.0
)";

    yaml node = yaml(contents);

    Database db = node;

    auto elements = db.elements();
    auto species = db.species();

    CHECK( elements.size() == 3 );
    CHECK( elements[0].symbol() == "A" );
    CHECK( elements[1].symbol() == "B" );
    CHECK( elements[2].symbol() == "C" );

    CHECK( species.size() == 5 );
    CHECK( species[0].name() == "A(aq)"      );
    CHECK( species[1].name() == "B(aq)"      );
    CHECK( species[2].name() == "A2B(aq)"    );
    CHECK( species[3].name() == "C(aq)"      );
    CHECK( species[4].name() == "A3B5C3(aq)" );

    CHECK( species[0].standardThermoModel().params().size() == 5 ); // G0, H0, V0, Cp0, Cv0
    CHECK( species[1].standardThermoModel().params().size() == 5 ); // G0, H0, V0, Cp0, Cv0
    CHECK( species[2].standardThermoModel().params().size() == 1 ); // lgKr
    CHECK( species[3].standardThermoModel().params().size() == 5 ); // G0, H0, V0, Cp0, Cv0
    CHECK( species[4].standardThermoModel().params().size() == 1 ); // lgKr

    CHECK( species[0].standardThermoModel().params()[0].value() == 1.0 );
    CHECK( species[1].standardThermoModel().params()[0].value() == 2.0 );
    CHECK( species[2].standardThermoModel().params()[0].value() == 5.0 );
    CHECK( species[3].standardThermoModel().params()[0].value() == 3.0 );
    CHECK( species[4].standardThermoModel().params()[0].value() == 7.0 );

    node = db; // convert back Database to yaml and check below for consistency

    CHECK( node["Elements"].size() == elements.size() );
    for(auto i = 0; i < elements.size(); ++i)
    {
        const auto enode = node["Elements"][i];
        CHECK( elements[i].symbol() == enode["Symbol"].as<String>() );
        CHECK( elements[i].molarMass() == enode["MolarMass"].as<double>() );
        CHECK( elements[i].name() == enode["Name"].as<String>() );
        if(enode["Tags"].IsDefined())
            CHECK( elements[i].tags() == enode["Tags"].as<Strings>() );
    }

    CHECK( node["Species"].size() == species.size() );
    for(auto i = 0; i < species.size(); ++i)
    {
        const auto snode = node["Species"][i];
        CHECK( species[i].name() == snode["Name"].as<String>() );
        CHECK( species[i].formula() == snode["Formula"].as<String>() );
        CHECK( species[i].substance() == snode["Substance"].as<String>() );
        CHECK( species[i].elements().repr() == snode["Elements"].as<String>() );
        CHECK( species[i].aggregateState() == snode["AggregateState"].as<AggregateState>() );
        if(snode["Tags"].IsDefined())
            CHECK( species[i].tags() == snode["Tags"].as<Strings>() );
    }
}

TEST_CASE("Testing YAML encoder/decoder for Element", "[Core.yaml]")
{
    Element H("H");

    yaml node;
    Element element;

    node = H;
    element = node;

    CHECK( element.symbol() == H.symbol() );
    CHECK( element.name() == H.name() );
    CHECK( element.molarMass() == H.molarMass() );
}

TEST_CASE("Testing YAML encoder/decoder for ElementList", "[Core.yaml]")
{
    yaml node;
    ElementList elementlist;

    // TODO: Implement YAML encoding/decoding test for ElementList.
}

TEST_CASE("Testing YAML encoder/decoder for ElementalComposition", "[Core.yaml]")
{
    yaml node;
    ElementalComposition elementalcomposition;

    // TODO: Implement YAML encoding/decoding test for ElementalComposition.
}

TEST_CASE("Testing YAML encoder/decoder for FormationReaction", "[Core.yaml]")
{
    yaml node;
    FormationReaction formationreaction;

    // TODO: Implement YAML encoding/decoding test for FormationReaction.
}

TEST_CASE("Testing YAML encoder/decoder for Param", "[Core.yaml]")
{
    yaml node;
    Param param;

    // TODO: Implement YAML encoding/decoding test for Param.
}

TEST_CASE("Testing YAML encoder/decoder for Params", "[Core.yaml]")
{
    yaml node;
    Params params;

    // TODO: Implement YAML encoding/decoding test for Params.
}

TEST_CASE("Testing YAML encoder/decoder for Phase", "[Core.yaml]")
{
    yaml node;
    Phase phase;

    // TODO: Implement YAML encoding/decoding test for Phase.
}

TEST_CASE("Testing YAML encoder/decoder for ReactionThermoModel", "[Core.yaml]")
{
    yaml node;
    ReactionThermoModel reactionthermomodel;

    // TODO: Implement YAML encoding/decoding test for ReactionThermoModel.
}

TEST_CASE("Testing YAML encoder/decoder for Species", "[Core.yaml]")
{
    yaml node;
    Species species;

    // TODO: Implement YAML encoding/decoding test for Species.
}

TEST_CASE("Testing YAML encoder/decoder for SpeciesList", "[Core.yaml]")
{
    yaml node;
    SpeciesList specieslist;

    // TODO: Implement YAML encoding/decoding test for SpeciesList.
}

TEST_CASE("Testing YAML encoder/decoder for StandardThermoModel", "[Core.yaml]")
{
    yaml node;
    StandardThermoModel standardthermomodel;

    // TODO: Implement YAML encoding/decoding test for StandardThermoModel.
}
