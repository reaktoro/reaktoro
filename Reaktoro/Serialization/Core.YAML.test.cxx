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
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardVolumeModelConstant.hpp>
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
      ReactionStandardThermoModel:
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
      ReactionStandardThermoModel:
        ConstLgK:
          lgKr: 7.0
)";

    yaml node = yaml::parse(contents);

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

    CHECK( species[0].standardThermoModel().params().size() == 6 ); // G0, H0, V0, Cp0, VT0, VP0
    CHECK( species[1].standardThermoModel().params().size() == 6 ); // G0, H0, V0, Cp0, VT0, VP0
    CHECK( species[2].standardThermoModel().params().size() == 1 ); // lgKr
    CHECK( species[3].standardThermoModel().params().size() == 6 ); // G0, H0, V0, Cp0, VT0, VP0
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
    ElementList elements = {
        Element("H"),
        Element("O"),
        Element("C"),
        Element("N")
    };

    yaml node = elements;

    for(auto i = 0; i < elements.size(); ++i)
    {
        CHECK( node[i]["Symbol"].IsDefined() );
        CHECK( node[i]["MolarMass"].IsDefined() );
        CHECK( node[i]["Name"].IsDefined() );

        CHECK( elements[i].symbol()    == node[i]["Symbol"].as<String>()    );
        CHECK( elements[i].molarMass() == node[i]["MolarMass"].as<double>() );
        CHECK( elements[i].name()      == node[i]["Name"].as<String>()      );

        if(node[i]["Tags"].IsDefined())
            CHECK( elements[i].tags() == node[i]["Tags"].as<Strings>()     );
    }

    ElementList elementlist = node;

    for(auto i = 0; i < elements.size(); ++i)
    {
        CHECK( elementlist[i].symbol() == elements[i].symbol() );
        CHECK( elementlist[i].molarMass() == elements[i].molarMass() );
        CHECK( elementlist[i].name() == elements[i].name() );
        CHECK( elementlist[i].tags() == elements[i].tags() );
    }
}

TEST_CASE("Testing YAML encoder/decoder for ElementalComposition", "[Core.yaml]")
{
    yaml node;
    ElementalComposition elements;

    elements = {{Element("H"), 2}, {Element("O"), 1}};
    node = elements;
    CHECK(node.as<String>() == "2:H 1:O");

    elements = {{Element("Ca"), 1}, {Element("C"), 1}, {Element("O"), 3}};
    node = elements;

    CHECK(node.as<String>() == "1:Ca 1:C 3:O");
}

TEST_CASE("Testing YAML encoder/decoder for FormationReaction", "[Core.yaml]")
{
    auto A = Species("Ca++").withStandardGibbsEnergy(0.0);
    auto B = Species("Mg++").withStandardGibbsEnergy(0.0);
    auto C = Species("CO3--").withStandardGibbsEnergy(0.0);

    auto reaction = FormationReaction()
        .withReactants({{A, 1}, {B, 1}, {C, 2}})
        .withEquilibriumConstant(1.0);

    yaml node = reaction;

    yaml expected = yaml::parse(R"(
Reactants: 1:Ca++ 1:Mg++ 2:CO3--
ReactionStandardThermoModel:
  ConstLgK:
    lgKr: 1
    Pr: 100000
StandardVolumeModel:
  Constant:
    V0: 0
)");

    CHECK( node.repr() == expected.repr() );
}

TEST_CASE("Testing YAML encoder/decoder for Param", "[Core.yaml]")
{
    Param param;
    yaml node;

    param = 1.0;
    node = param;

    CHECK( node.as<double>() == 1.0 );

    node = 10.0;
    param = node.as<double>();

    CHECK( param.value() == 10.0 );
}

TEST_CASE("Testing YAML encoder/decoder for Vec<Param>", "[Core.yaml]")
{
    yaml node = Vec<double>{1.0, 2.0, 3.0};

    Vec<Param> params = node;

    CHECK( params[0].value() == 1.0 );
    CHECK( params[1].value() == 2.0 );
    CHECK( params[2].value() == 3.0 );

    params[0] = 10.0;
    params[1] = 20.0;
    params[2] = 30.0;

    node = params;

    CHECK( node[0].as<double>() == 10.0 );
    CHECK( node[1].as<double>() == 20.0 );
    CHECK( node[2].as<double>() == 30.0 );
}

TEST_CASE("Testing YAML encoder/decoder for Phase", "[Core.yaml]")
{
    yaml node;
    Phase phase;

    // TODO: Implement YAML encoding/decoding test for Phase.
}

TEST_CASE("Testing YAML encoder/decoder for ReactionStandardThermoModel", "[Core.yaml]")
{
    yaml node, expected;

    Param lgKr = 1.0;
    Param dHr  = 2.0;
    Param Tr   = 3.0;
    Param Pr   = 4.0;

    node = ReactionStandardThermoModelVantHoff({lgKr, dHr, Tr, Pr});

    expected = yaml::parse(R"(
        VantHoff:
          lgKr: 1
          dHr: 2
          Tr: 3
          Pr: 4
    )");

    CHECK( node.repr() == expected.repr() );

    auto model = ReactionStandardThermoModelYAML(node);

    CHECK( model.serialize().repr() == expected.repr() );
}

TEST_CASE("Testing YAML encoder/decoder for Species", "[Core.yaml]")
{
    Species species;
    yaml node, expected;

    WHEN("using constructor Species(formula)")
    {
        node = Species("CaCO3(aq)")
            .withStandardGibbsEnergy(10.0);

        expected = yaml::parse(R"(
            Name: CaCO3(aq)
            Formula: CaCO3
            Substance: CaCO3
            Elements: 1:Ca 1:C 3:O
            AggregateState: Aqueous
            StandardThermoModel:
              Constant:
                G0: 10
                H0: 0
                V0: 0
                VT0: 0
                VP0: 0
                Cp0: 0
        )");

        CHECK( node.repr() == expected.repr() );
    }

    WHEN("using constructor Species(formula) with HKF standard thermodynamic model and charged species")
    {
        StandardThermoModelParamsHKF params;
        params.Gf = 1.0;
        params.Hf = 2.0;
        params.Sr = 3.0;
        params.a1 = 4.0;
        params.a2 = 5.0;
        params.a3 = 6.0;
        params.a4 = 7.0;
        params.c1 = 8.0;
        params.c2 = 9.0;
        params.wref = 10.0;
        params.charge = 11.0;
        params.Tmax = 12.0;

        node = Species("CO3--(aq)")
            .withSubstance("CARBONATE")
            .withStandardThermoModel(StandardThermoModelHKF(params));

        expected = yaml::parse(R"(
            Name: CO3--(aq)
            Formula: CO3--
            Substance: CARBONATE
            Elements: 1:C 3:O
            Charge: -2
            AggregateState: Aqueous
            StandardThermoModel:
              HKF:
                Gf: 1
                Hf: 2
                Sr: 3
                a1: 4
                a2: 5
                a3: 6
                a4: 7
                c1: 8
                c2: 9
                wref: 10
                charge: 11
                Tmax: 12
        )");

        CHECK( node.repr() == expected.repr() );
    }

    WHEN("using constructor Species(formula) with FormationReaction")
    {
        auto A = Species("Ca++(aq)").withStandardGibbsEnergy(0.0);
        auto B = Species("CO3--(aq)").withStandardGibbsEnergy(0.0);

        Param lgKr = 1.0;
        Param dHr  = 2.0;
        Param Tr   = 3.0;
        Param Pr   = 4.0;
        Param V0   = 5.0;

        auto reaction = FormationReaction()
            .withReactants({{A, 1}, {B, 1}})
            .withReactionStandardThermoModel(ReactionStandardThermoModelVantHoff({lgKr, dHr, Tr, Pr}))
            .withProductStandardVolumeModel(StandardVolumeModelConstant({V0}));

        node = Species("CaCO3(s)").withFormationReaction(reaction);

        expected = yaml::parse(R"(
            Name: CaCO3(s)
            Formula: CaCO3
            Substance: CaCO3
            Elements: 1:Ca 1:C 3:O
            AggregateState: Solid
            FormationReaction:
              Reactants: 1:Ca++(aq) 1:CO3--(aq)
              ReactionStandardThermoModel:
                VantHoff:
                  lgKr: 1
                  dHr: 2
                  Tr: 3
                  Pr: 4
              StandardVolumeModel:
                Constant:
                  V0: 5
        )");

        CHECK( node.repr() == expected.repr() );
    }
}

TEST_CASE("Testing YAML encoder/decoder for SpeciesList", "[Core.yaml]")
{
    yaml node = SpeciesList({
        Species("Ca++(aq)").withStandardGibbsEnergy(0.0),
        Species("CO3--(aq)").withStandardGibbsEnergy(0.0),
        Species("CaCO3(aq)").withStandardGibbsEnergy(0.0)
    });

    yaml expected = yaml::parse(R"(
        - Name: Ca++(aq)
          Formula: Ca++
          Substance: Ca++
          Elements: 1:Ca
          Charge: 2
          AggregateState: Aqueous
          StandardThermoModel:
            Constant:
              G0: 0
              H0: 0
              V0: 0
              VT0: 0
              VP0: 0
              Cp0: 0
        - Name: CO3--(aq)
          Formula: CO3--
          Substance: CO3--
          Elements: 1:C 3:O
          Charge: -2
          AggregateState: Aqueous
          StandardThermoModel:
            Constant:
              G0: 0
              H0: 0
              V0: 0
              VT0: 0
              VP0: 0
              Cp0: 0
        - Name: CaCO3(aq)
          Formula: CaCO3
          Substance: CaCO3
          Elements: 1:Ca 1:C 3:O
          AggregateState: Aqueous
          StandardThermoModel:
            Constant:
              G0: 0
              H0: 0
              V0: 0
              VT0: 0
              VP0: 0
              Cp0: 0
    )");

    CHECK( node.repr() == expected.repr() );
}

TEST_CASE("Testing YAML encoder/decoder for StandardThermoModel", "[Core.yaml]")
{
    StandardThermoModelParamsMaierKelley params;
    params.Gf = 1.0;
    params.Hf = 2.0;
    params.Sr = 3.0;
    params.Vr = 4.0;
    params.a = 5.0;
    params.b = 6.0;
    params.c = 7.0;
    params.Tmax = 8.0;

    yaml node = StandardThermoModelMaierKelley(params);

    yaml expected = yaml::parse(R"(
        MaierKelley:
          Gf: 1
          Hf: 2
          Sr: 3
          Vr: 4
          a: 5
          b: 6
          c: 7
          Tmax: 8
    )");

    CHECK( node.repr() == expected.repr() );

    auto model = StandardThermoModelYAML(node);

    CHECK( model.serialize().repr() == expected.repr() );
}
