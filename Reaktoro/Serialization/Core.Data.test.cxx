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
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Serialization/Common.Data.hpp>
#include <Reaktoro/Serialization/Core.Data.hpp>
#include <Reaktoro/Core/Data.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelYAML.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardVolumeModelConstant.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Data encoder/decoder for AggregateState", "[Core.Data]")
{
    Data data;
    AggregateState aggstate;

    data = AggregateState::Aqueous;
    CHECK( data.asString() == "Aqueous" );
    aggstate = data.as<AggregateState>();
    CHECK( aggstate == AggregateState::Aqueous );

    data = AggregateState::Gas;
    CHECK( data.asString() == "Gas" );
    aggstate = data.as<AggregateState>();
    CHECK( aggstate == AggregateState::Gas );

    data = AggregateState::Solid;
    CHECK( data.asString() == "Solid" );
    aggstate = data.as<AggregateState>();
    CHECK( aggstate == AggregateState::Solid );

    data = AggregateState::Undefined;
    CHECK( data.asString() == "Undefined" );
    aggstate = data.as<AggregateState>();
    CHECK( aggstate == AggregateState::Undefined );
}

TEST_CASE("Testing Data encoder/decoder for ChemicalFormula", "[Core.Data]")
{
    Data data;
    ChemicalFormula formula;

    data = ChemicalFormula("H2O");
    CHECK( data.asString() == "H2O" );
    formula = data.as<ChemicalFormula>();
    CHECK( formula.equivalent("H2O") );

    data = ChemicalFormula("Ca++");
    CHECK( data.asString() == "Ca++" );
    formula = data.as<ChemicalFormula>();
    CHECK( formula.equivalent("Ca++") );
}

TEST_CASE("Testing Data encoder/decoder for ChemicalSystem", "[Core.Data]")
{
    Data data;
    ChemicalSystem chemicalsystem;

    // TODO: Implement tests for Data encoding/decoding of ChemicalSystem objects.
}

TEST_CASE("Testing Data encoder/decoder for Database", "[Core.Data]")
{
    auto contents = R"(
Elements:
  A:
    Symbol: A
    MolarMass: 1.0
  B:
    Symbol: B
    MolarMass: 2.0
  C:
    Symbol: C
    MolarMass: 3.0
Species:
  A2B(aq):
    Name: A2B(aq)
    Formula: A2B
    Elements: 2:A 1:B
    AggregateState: Aqueous
    FormationReaction:
      Reactants: 2:A(aq) 1:B(aq)
      ReactionStandardThermoModel:
        ConstLgK:
          lgKr: 5.0
  A(aq):
    Name: A(aq)
    Formula: A
    Elements: 1:A
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 1.0 }
  B(aq):
    Name: B(aq)
    Formula: B
    Elements: 1:B
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 2.0 }
  C(aq):
    Name: C(aq)
    Formula: C
    Elements: 1:C
    AggregateState: Aqueous
    StandardThermoModel:
      Constant: { G0: 3.0 }
  A3B5C3(aq):
    Name: A3B5C3(aq)
    Formula: A3B5C3
    Elements: 3:A 5:B 3:C
    AggregateState: Aqueous
    FormationReaction:
      Reactants: 1:A2B(aq) 1:A(aq) 4:B(aq) 3:C(aq)
      ReactionStandardThermoModel:
        ConstLgK:
          lgKr: 7.0
)";

    Data data = Data::fromYaml(contents);

    Database db = data.as<Database>();

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

    data = db; // convert back Database to Data and check below for consistency

    CHECK( data["Elements"].asList().size() == elements.size() );
    for(auto i = 0; i < elements.size(); ++i)
    {
        const auto enode = data["Elements"][elements[i].symbol()];
        CHECK( elements[i].symbol() == enode["Symbol"].asString() );
        CHECK( elements[i].molarMass() == enode["MolarMass"].asFloat() );
        CHECK( elements[i].name() == enode["Name"].asString() );
        if(enode.exists("Tags"))
            CHECK( elements[i].tags() == enode["Tags"].as<Strings>() );
    }

    CHECK( data["Species"].asList().size() == species.size() );
    for(auto i = 0; i < species.size(); ++i)
    {
        const auto snode = data["Species"][species[i].name()];
        CHECK( species[i].name() == snode["Name"].asString() );
        CHECK( species[i].formula() == snode["Formula"].asString() );
        CHECK( species[i].substance() == snode["Substance"].asString() );
        CHECK( species[i].elements().repr() == snode["Elements"].asString() );
        CHECK( species[i].aggregateState() == snode["AggregateState"].as<AggregateState>() );
        if(snode.exists("Tags"))
            CHECK( species[i].tags() == snode["Tags"].as<Strings>() );
    }
}

TEST_CASE("Testing Data encoder/decoder for Element", "[Core.Data]")
{
    Element H("H");

    Data data;
    Element element;

    data = H;
    element = data.as<Element>();

    CHECK( element.symbol() == H.symbol() );
    CHECK( element.name() == H.name() );
    CHECK( element.molarMass() == H.molarMass() );
}

TEST_CASE("Testing Data encoder/decoder for ElementList", "[Core.Data]")
{
    ElementList elements = {
        Element("H"),
        Element("O"),
        Element("C"),
        Element("N")
    };

    Data data = elements;

    for(auto i = 0; i < elements.size(); ++i)
    {
        CHECK( data[i].exists("Symbol") );
        CHECK( data[i].exists("MolarMass") );
        CHECK( data[i].exists("Name") );

        CHECK( elements[i].symbol()    == data[i]["Symbol"].asString()   );
        CHECK( elements[i].molarMass() == data[i]["MolarMass"].asFloat() );
        CHECK( elements[i].name()      == data[i]["Name"].asString()     );

        if(data[i].exists("Tags"))
            CHECK( elements[i].tags() == data[i]["Tags"].as<Strings>()     );
    }

    ElementList elementlist = data.as<ElementList>();

    for(auto i = 0; i < elements.size(); ++i)
    {
        CHECK( elementlist[i].symbol() == elements[i].symbol() );
        CHECK( elementlist[i].molarMass() == elements[i].molarMass() );
        CHECK( elementlist[i].name() == elements[i].name() );
        CHECK( elementlist[i].tags() == elements[i].tags() );
    }
}

TEST_CASE("Testing Data encoder/decoder for ElementalComposition", "[Core.Data]")
{
    Data data;
    ElementalComposition elements;

    elements = {{Element("H"), 2}, {Element("O"), 1}};
    data = elements;
    CHECK(data.asString() == "2:H 1:O");

    elements = {{Element("Ca"), 1}, {Element("C"), 1}, {Element("O"), 3}};
    data = elements;

    CHECK(data.asString() == "1:Ca 1:C 3:O");
}

TEST_CASE("Testing Data encoder/decoder for FormationReaction", "[Core.Data]")
{
    auto A = Species("Ca++").withStandardGibbsEnergy(0.0);
    auto B = Species("Mg++").withStandardGibbsEnergy(0.0);
    auto C = Species("CO3--").withStandardGibbsEnergy(0.0);

    auto reaction = FormationReaction()
        .withReactants({{A, 1}, {B, 1}, {C, 2}})
        .withEquilibriumConstant(1.0);

    Data data = reaction;

    Data expected = Data::fromYaml(R"(
Reactants: 1:Ca++ 1:Mg++ 2:CO3--
ReactionStandardThermoModel:
  ConstLgK:
    lgKr: 1
    Pr: 100000
StandardVolumeModel:
  Constant:
    V0: 0
)");

    CHECK( data.dumpYaml() == expected.dumpYaml() );
}

TEST_CASE("Testing Data encoder/decoder for Param", "[Core.Data]")
{
    Param param;
    Data data;

    param = 1.0;
    data = param;

    CHECK( data.asFloat() == 1.0 );

    data = 10.0;
    param = data.asFloat();

    CHECK( param.value() == 10.0 );
}

TEST_CASE("Testing Data encoder/decoder for Vec<Param>", "[Core.Data]")
{
    Data data = Vec<double>{1.0, 2.0, 3.0};

    Vec<Param> params = data.as<Vec<Param>>();

    CHECK( params[0].value() == 1.0 );
    CHECK( params[1].value() == 2.0 );
    CHECK( params[2].value() == 3.0 );

    params[0] = 10.0;
    params[1] = 20.0;
    params[2] = 30.0;

    data = params;

    CHECK( data[0].asFloat() == 10.0 );
    CHECK( data[1].asFloat() == 20.0 );
    CHECK( data[2].asFloat() == 30.0 );
}

TEST_CASE("Testing Data encoder/decoder for Phase", "[Core.Data]")
{
    Data data;
    Phase phase;

    // TODO: Implement Data encoding/decoding test for Phase.
}

TEST_CASE("Testing Data encoder/decoder for ReactionStandardThermoModel", "[Core.Data]")
{
    Data data, expected;

    Param lgKr = 1.0;
    Param dHr  = 2.0;
    Param Tr   = 3.0;
    Param Pr   = 4.0;

    data = ReactionStandardThermoModelVantHoff({lgKr, dHr, Tr, Pr});

    expected = Data::fromYaml(R"(
        VantHoff:
          lgKr: 1
          dHr: 2
          Tr: 3
          Pr: 4
    )");

    CHECK( data.dumpYaml() == expected.dumpYaml() );

    // auto model = ReactionStandardThermoModelYAML(data);

    // CHECK( model.serialize().dumpYaml() == expected.dumpYaml() );
}

TEST_CASE("Testing Data encoder/decoder for Species", "[Core.Data]")
{
    Species species;
    Data data, expected;

    WHEN("using constructor Species(formula)")
    {
        data = Species("CaCO3(aq)")
            .withStandardGibbsEnergy(10.0);

        expected = Data::fromYaml(R"(
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

        CHECK( data.dumpYaml() == expected.dumpYaml() );
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

        data = Species("CO3--(aq)")
            .withSubstance("CARBONATE")
            .withStandardThermoModel(StandardThermoModelHKF(params));

        expected = Data::fromYaml(R"(
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

        CHECK( data.dumpYaml() == expected.dumpYaml() );
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

        data = Species("CaCO3(s)").withFormationReaction(reaction);

        expected = Data::fromYaml(R"(
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

        CHECK( data.dumpYaml() == expected.dumpYaml() );
    }
}

TEST_CASE("Testing Data encoder/decoder for SpeciesList", "[Core.Data]")
{
    Data data = SpeciesList({
        Species("Ca++(aq)").withStandardGibbsEnergy(0.0),
        Species("CO3--(aq)").withStandardGibbsEnergy(0.0),
        Species("CaCO3(aq)").withStandardGibbsEnergy(0.0)
    });

    Data expected = Data::fromYaml(R"(
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

    CHECK( data.dumpYaml() == expected.dumpYaml() );
}

TEST_CASE("Testing Data encoder/decoder for StandardThermoModel", "[Core.Data]")
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

    Data data = StandardThermoModelMaierKelley(params);

    Data expected = Data::fromYaml(R"(
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

    CHECK( data.dumpYaml() == expected.dumpYaml() );

    // auto model = StandardThermoModelYAML(data);

    // CHECK( model.serialize().dumpYaml() == expected.dumpYaml() );
}
