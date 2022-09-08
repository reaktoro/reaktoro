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

// YAML/JSON includes
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>

// Reaktoro includes
#include <Reaktoro/Core/Data.hpp>
using namespace Reaktoro;

const auto yaml_testing_string = R"#(
Species:
- Name: Almandine
  Formula: Fe3Al2Si3O12
  Elements: 3:Fe 2:Al 3:Si 12:O
  AggregateState: Solid
  StandardThermoModel:
    HollandPowell:
      Gf: -4937500.0
      Hf: -5260650.0
      Sr: 342.0
      Vr: 0.00011525
      a: 677.3
      b: 0.0
      c: -3772700.0
      d: -5044.0
      alpha0: 2.12e-05
      kappa0: 190000000000.0
      kappa0p: 2.98
      kappa0pp: -1.6e-11
      numatoms: 20.0
      Tmax: 9999.0
- Name: Andradite
  Formula: Ca3Fe2Si3O12
  Elements: 3:Ca 2:Fe 3:Si 12:O
  AggregateState: Solid
  StandardThermoModel:
    HollandPowell:
      Gf: -5426110.0
      Hf: -5769080.0
      Sr: 316.4
      Vr: 0.00013204
      a: 638.6
      b: 0.0
      c: -4955100.0
      d: -3989.2
      alpha0: 2.86e-05
      kappa0: 158800000000.0
      kappa0p: 5.68
      kappa0pp: -3.6e-11
      numatoms: 20.0
      Tmax: 9999.0
Extra:
    SomeBoolean: true
    AnotherBoolean: false
    AnInteger: 123
    SomeStrings: [Hello, Hallo]
)#";

const auto json_testing_string = R"#(
{
  "Species": [
    {
      "Name": "Almandine",
      "Formula": "Fe3Al2Si3O12",
      "Elements": "3:Fe 2:Al 3:Si 12:O",
      "AggregateState": "Solid",
      "StandardThermoModel": {
        "HollandPowell": {
          "Gf": -4937500.0,
          "Hf": -5260650.0,
          "Sr": 342.0,
          "Vr": 0.00011525,
          "a": 677.3,
          "b": 0.0,
          "c": -3772700.0,
          "d": -5044.0,
          "alpha0": 2.12e-05,
          "kappa0": 190000000000.0,
          "kappa0p": 2.98,
          "kappa0pp": -1.6e-11,
          "numatoms": 20.0,
          "Tmax": 9999.0
        }
      }
    },
    {
      "Name": "Andradite",
      "Formula": "Ca3Fe2Si3O12",
      "Elements": "3:Ca 2:Fe 3:Si 12:O",
      "AggregateState": "Solid",
      "StandardThermoModel": {
        "HollandPowell": {
          "Gf": -5426110.0,
          "Hf": -5769080.0,
          "Sr": 316.4,
          "Vr": 0.00013204,
          "a": 638.6,
          "b": 0.0,
          "c": -4955100.0,
          "d": -3989.2,
          "alpha0": 2.86e-05,
          "kappa0": 158800000000.0,
          "kappa0p": 5.68,
          "kappa0pp": -3.6e-11,
          "numatoms": 20.0,
          "Tmax": 9999.0
        }
      }
    }
  ],
  "Extra": {
    "SomeBoolean": true,
    "AnotherBoolean": false,
    "AnInteger": 123,
    "SomeStrings": ["Hello", "Hallo"]
  }
}
)#";

// Used to test encoding/decoding of custom types to/from Data objects
struct ComplicatedStuff
{
    bool boolval;               // required
    double doubleval;           // required
    Param paramval;             // required
    String stringval;           // required
    Strings vecval;             // optional
    Map<String, double> mapval; // optional
};

REAKTORO_DATA_ENCODE_DECLARE(ComplicatedStuff);
REAKTORO_DATA_DECODE_DECLARE(ComplicatedStuff);

REAKTORO_DATA_ENCODE_DEFINE(ComplicatedStuff)
{
    data["boolval"] = obj.boolval;
    data["doubleval"] = obj.doubleval;
    data["paramval"] = obj.paramval;
    data["stringval"] = obj.stringval;
    data["vecval"] = obj.vecval;
    data["mapval"] = obj.mapval;
}

REAKTORO_DATA_DECODE_DEFINE(ComplicatedStuff)
{
    data.required("boolval").to(obj.boolval);
    data.required("doubleval").to(obj.doubleval);
    data.required("paramval").to(obj.paramval);
    data.required("stringval").to(obj.stringval);
    data.optional("vecval").to(obj.vecval);
    data.optional("mapval").to(obj.mapval);
}

TEST_CASE("Testing Data class", "[Data]")
{
    SECTION("Checking Data::add(value) for the creation of lists")
    {
        Data data;

        // Make a list by using Data::add(value)
        data.add(true);
        data.add(false);
        data.add("ABC");
        data.add(String("CDE"));
        data.add(1);
        data.add(-2);
        data.add(std::int64_t(-3));
        data.add(std::uint64_t(4));
        data.add(std::float_t(5.0));
        data.add(std::double_t(6.0));
        data.add(real(7.0));
        data.add(Param(8.0));
        data.add(Vec<Data>{data[0], data[1]});
        data.add(Dict<String, Data>{{"U", data[2]}, {"V", data[3]}});
        data.add(nullptr);

        // Check the types of the Data objects created
        CHECK( data[0].isBoolean() );
        CHECK( data[1].isBoolean() );
        CHECK( data[2].isString()  );
        CHECK( data[3].isString()  );
        CHECK( data[4].isParam()   );
        CHECK( data[5].isParam()   );
        CHECK( data[6].isParam()   );
        CHECK( data[7].isParam()   );
        CHECK( data[8].isParam()   );
        CHECK( data[9].isParam()   );
        CHECK( data[10].isParam()  );
        CHECK( data[11].isParam()  );
        CHECK( data[12].isList()   );
        CHECK( data[13].isDict()   );
        CHECK( data[14].isNull()   );

        // Check the values of the Data objects created via Data::add(value)
        CHECK( data[0].asBoolean()      ==  true    );
        CHECK( data[1].asBoolean()      ==  false   );
        CHECK( data[2].asString()       ==  "ABC"   );
        CHECK( data[3].asString()       ==  "CDE"   );
        CHECK( data[4].asParam()        ==  1.0     );
        CHECK( data[5].asParam()        == -2.0     );
        CHECK( data[6].asParam()        == -3.0     );
        CHECK( data[7].asParam()        ==  4.0     );
        CHECK( data[8].asParam()        ==  5.0     );
        CHECK( data[9].asParam()        ==  6.0     );
        CHECK( data[10].asParam()       ==  7.0     );
        CHECK( data[11].asParam()       ==  8.0     );
        CHECK( data[12][0].asBoolean()  ==  true    );
        CHECK( data[12][1].asBoolean()  ==  false   );
        CHECK( data[13]["U"].asString() ==  "ABC"   );
        CHECK( data[13]["V"].asString() ==  "CDE"   );
        CHECK( data[14].asNull()        ==  nullptr );

        // Modify data[0], data[5] and data[9]
        data[0].assign(1.0);
        data[5].assign("Hello");
        data[9].reset();

        CHECK( data[0].isParam()  );
        CHECK( data[5].isString() );
        CHECK( data[9].isNull()   );

        CHECK( data[0].asParam()  == 1.0     );
        CHECK( data[5].asString() == "Hello" );
        CHECK( data[9].asNull()   == nullptr );
    }

    SECTION("Checking Data::add(key, value) for the creation of dictionaries")
    {
        Data data;

        // Make a dictionary by using Data::add(key, value)
        data.add("K0", true);
        data.add("K1", false);
        data.add("K2", "ABC");
        data.add("K3", String("CDE"));
        data.add("K4", 1);
        data.add("K5", -2);
        data.add("K6", std::int64_t(-3));
        data.add("K7", std::uint64_t(4));
        data.add("K8", std::float_t(5.0));
        data.add("K9", std::double_t(6.0));
        data.add("K10", real(7.0));
        data.add("K11", Param(8.0));
        data.add("K12", Vec<Data>{data["K0"], data["K1"]});
        data.add("K13", Dict<String, Data>{{"U", data["K2"]}, {"V", data["K3"]}});
        data.add("K14", nullptr);

        // Check the types of the Data objects created
        CHECK( data["K0"].isBoolean() );
        CHECK( data["K1"].isBoolean() );
        CHECK( data["K2"].isString()  );
        CHECK( data["K3"].isString()  );
        CHECK( data["K4"].isParam()   );
        CHECK( data["K5"].isParam()   );
        CHECK( data["K6"].isParam()   );
        CHECK( data["K7"].isParam()   );
        CHECK( data["K8"].isParam()   );
        CHECK( data["K9"].isParam()   );
        CHECK( data["K10"].isParam()  );
        CHECK( data["K11"].isParam()  );
        CHECK( data["K12"].isList()   );
        CHECK( data["K13"].isDict()   );
        CHECK( data["K14"].isNull()   );

        // Check the values of the Data objects created via Data::add(key, value)
        CHECK( data["K0"].asBoolean()      ==  true    );
        CHECK( data["K1"].asBoolean()      ==  false   );
        CHECK( data["K2"].asString()       ==  "ABC"   );
        CHECK( data["K3"].asString()       ==  "CDE"   );
        CHECK( data["K4"].asParam()        ==  1.0     );
        CHECK( data["K5"].asParam()        == -2.0     );
        CHECK( data["K6"].asParam()        == -3.0     );
        CHECK( data["K7"].asParam()        ==  4.0     );
        CHECK( data["K8"].asParam()        ==  5.0     );
        CHECK( data["K9"].asParam()        ==  6.0     );
        CHECK( data["K10"].asParam()       ==  7.0     );
        CHECK( data["K11"].asParam()       ==  8.0     );
        CHECK( data["K12"][0].asBoolean()  ==  true    );
        CHECK( data["K12"][1].asBoolean()  ==  false   );
        CHECK( data["K13"]["U"].asString() ==  "ABC"   );
        CHECK( data["K13"]["V"].asString() ==  "CDE"   );
        CHECK( data["K14"].asNull()        ==  nullptr );

        // Modify data[0], data[5] and data[9]
        data["K0"].assign(1.0);
        data["K5"].assign("Hello");
        data["K9"].reset();

        CHECK( data["K0"].isParam()  );
        CHECK( data["K5"].isString() );
        CHECK( data["K9"].isNull()   );

        CHECK( data["K0"].asParam()  == 1.0     );
        CHECK( data["K5"].asString() == "Hello" );
        CHECK( data["K9"].asNull()   == nullptr );
    }

    SECTION("Checking the construction of nested Data objects using operator[]")
    {
        Data foo;
        foo["A"] = 1.0;
        foo["B"] = 2;

        Data bar;
        bar["C"] = true;
        bar["D"] = "X";
        bar["E"] = Param(7.0);

        Data vec = Vec<Data>{ Data(false), Data("Y"), Data(9.0) };

        Data data;
        data["Foo"] = foo;
        data["Bar"] = bar;
        data["Vec"] = vec;

        SECTION("Check if data inserted is correct")
        {
            CHECK( data.exists("Foo") == true  );
            CHECK( data.exists("Bar") == true  );
            CHECK( data.exists("Vec") == true  );
            CHECK( data.exists("Joe") == false );

            CHECK( data["Foo"].exists("A") == true  );
            CHECK( data["Foo"].exists("B") == true  );
            CHECK( data["Foo"].exists("Z") == false );

            CHECK( data["Bar"].exists("C") == true  );
            CHECK( data["Bar"].exists("D") == true  );
            CHECK( data["Bar"].exists("E") == true  );
            CHECK( data["Bar"].exists("Z") == false );

            CHECK( data["Vec"].isList() );
            CHECK( data["Vec"].asList().size() == 3 );

            CHECK( data["Foo"]["A"].asFloat()   == 1.0  );
            CHECK( data["Foo"]["B"].asInteger() == 2    );

            CHECK( data["Bar"]["C"].asBoolean() == true );
            CHECK( data["Bar"]["D"].asString()  == "X"  );
            CHECK( data["Bar"]["E"].asParam()   == 7.0  );

            CHECK( data["Vec"][0].asBoolean() == false );
            CHECK( data["Vec"][1].asString()  == "Y"   );
            CHECK( data["Vec"][2].asParam()   == 9.0   );
        }

        SECTION("Check exceptions on bad use of accessor methods")
        {
            CHECK_THROWS( data.at("Joe") );
            CHECK_THROWS( data.required("Joe") );

            CHECK_NOTHROW( data.optional("Joe") );
            CHECK_NOTHROW( data["Roe"] );

            CHECK( data.exists("Roe") ); // Roe introduced before with operator[]
        }
    }

    SECTION("Check conversion from integer and floating point types to Param objects and vice-versa")
    {
        const auto num0 = Data(1);
        const auto num1 = Data(-2);
        const auto num2 = Data(std::int64_t(-3));
        const auto num3 = Data(std::uint64_t(4));
        const auto num4 = Data(std::float_t(5.2));
        const auto num5 = Data(std::double_t(6.4));
        const auto num6 = Data(real(7.8));
        const auto num7 = Data(Param(8.3));
        const auto num8 = Data(true);
        const auto num9 = Data('a');

        CHECK( num0.isParam() == true  );
        CHECK( num1.isParam() == true  );
        CHECK( num2.isParam() == true  );
        CHECK( num3.isParam() == true  );
        CHECK( num4.isParam() == true  );
        CHECK( num5.isParam() == true  );
        CHECK( num6.isParam() == true  );
        CHECK( num7.isParam() == true  );
        CHECK( num8.isParam() == false ); // num8 should be boolean
        CHECK( num9.isParam() == false ); // num9 should be string

        CHECK( num0.asInteger() ==  1 );
        CHECK( num1.asInteger() == -2 );
        CHECK( num2.asInteger() == -3 );
        CHECK( num3.asInteger() ==  4 );
        CHECK( num4.asInteger() ==  5 );
        CHECK( num5.asInteger() ==  6 );
        CHECK( num6.asInteger() ==  7 );
        CHECK( num7.asInteger() ==  8 );

        CHECK( num0.asFloat() ==  Approx( 1.0) );
        CHECK( num1.asFloat() ==  Approx(-2.0) );
        CHECK( num2.asFloat() ==  Approx(-3.0) );
        CHECK( num3.asFloat() ==  Approx( 4.0) );
        CHECK( num4.asFloat() ==  Approx( 5.2) );
        CHECK( num5.asFloat() ==  Approx( 6.4) );
        CHECK( num6.asFloat() ==  Approx( 7.8) );
        CHECK( num7.asFloat() ==  Approx( 8.3) );

        CHECK( num0.asParam().value().val() ==  Approx( 1.0) );
        CHECK( num1.asParam().value().val() ==  Approx(-2.0) );
        CHECK( num2.asParam().value().val() ==  Approx(-3.0) );
        CHECK( num3.asParam().value().val() ==  Approx( 4.0) );
        CHECK( num4.asParam().value().val() ==  Approx( 5.2) );
        CHECK( num5.asParam().value().val() ==  Approx( 6.4) );
        CHECK( num6.asParam().value().val() ==  Approx( 7.8) );
        CHECK( num7.asParam().value().val() ==  Approx( 8.3) );

        CHECK_NOTHROW( num8.asBoolean() );
        CHECK_NOTHROW( num9.asString()  );

        CHECK_THROWS( num8.asString()   );
        CHECK_THROWS( num8.asInteger()  );
        CHECK_THROWS( num8.asFloat()    );
        CHECK_THROWS( num8.asReal()     );
        CHECK_THROWS( num8.asParam()    );
        CHECK_THROWS( num8.asDict()     );
        CHECK_THROWS( num8.asList()     );
        CHECK_THROWS( num8.asNull()     );

        CHECK_THROWS( num9.asBoolean()  );
        CHECK_THROWS( num9.asInteger()  );
        CHECK_THROWS( num9.asFloat()    );
        CHECK_THROWS( num9.asReal()     );
        CHECK_THROWS( num9.asParam()    );
        CHECK_THROWS( num9.asDict()     );
        CHECK_THROWS( num9.asList()     );
        CHECK_THROWS( num9.asNull()     );
    }

    SECTION("Checking the construction of Data objects using YAML and JSON formatted strings")
    {
        const Data data = GENERATE(
            Data::parseYaml(yaml_testing_string),
            Data::parseJson(json_testing_string)
        );

        CHECK( data.isDict() );

        CHECK( data["Species"].isList() );
        CHECK( data["Species"].asList().size() == 2 );

        const Data species0 = data["Species"][0];
        const Data species1 = data["Species"][1];

        CHECK( species0["Name"].asString()           == "Almandine"           );
        CHECK( species0["Formula"].asString()        == "Fe3Al2Si3O12"        );
        CHECK( species0["Elements"].asString()       == "3:Fe 2:Al 3:Si 12:O" );
        CHECK( species0["AggregateState"].asString() == "Solid"               );

        CHECK( species1["Name"].asString()           == "Andradite"           );
        CHECK( species1["Formula"].asString()        == "Ca3Fe2Si3O12"        );
        CHECK( species1["Elements"].asString()       == "3:Ca 2:Fe 3:Si 12:O" );
        CHECK( species1["AggregateState"].asString() == "Solid"               );

        const Data thermo0 = species0["StandardThermoModel"]["HollandPowell"];
        const Data thermo1 = species1["StandardThermoModel"]["HollandPowell"];

        CHECK( thermo0["Gf"].asFloat()       == -4937500.0      );
        CHECK( thermo0["Hf"].asFloat()       == -5260650.0      );
        CHECK( thermo0["Sr"].asFloat()       ==  342.0          );
        CHECK( thermo0["Vr"].asFloat()       ==  0.00011525     );
        CHECK( thermo0["a"].asFloat()        ==  677.3          );
        CHECK( thermo0["b"].asFloat()        ==  0.0            );
        CHECK( thermo0["c"].asFloat()        == -3772700.0      );
        CHECK( thermo0["d"].asFloat()        == -5044.0         );
        CHECK( thermo0["alpha0"].asFloat()   ==  2.12e-05       );
        CHECK( thermo0["kappa0"].asFloat()   ==  190000000000.0 );
        CHECK( thermo0["kappa0p"].asFloat()  ==  2.98           );
        CHECK( thermo0["kappa0pp"].asFloat() == -1.6e-11        );
        CHECK( thermo0["numatoms"].asFloat() ==  20.0           );
        CHECK( thermo0["Tmax"].asFloat()     ==  9999.0         );

        CHECK( thermo1["Gf"].asFloat()       == -5426110.0      );
        CHECK( thermo1["Hf"].asFloat()       == -5769080.0      );
        CHECK( thermo1["Sr"].asFloat()       ==  316.4          );
        CHECK( thermo1["Vr"].asFloat()       ==  0.00013204     );
        CHECK( thermo1["a"].asFloat()        ==  638.6          );
        CHECK( thermo1["b"].asFloat()        ==  0.0            );
        CHECK( thermo1["c"].asFloat()        == -4955100.0      );
        CHECK( thermo1["d"].asFloat()        == -3989.2         );
        CHECK( thermo1["alpha0"].asFloat()   ==  2.86e-05       );
        CHECK( thermo1["kappa0"].asFloat()   ==  158800000000.0 );
        CHECK( thermo1["kappa0p"].asFloat()  ==  5.68           );
        CHECK( thermo1["kappa0pp"].asFloat() == -3.6e-11        );
        CHECK( thermo1["numatoms"].asFloat() ==  20.0           );
        CHECK( thermo1["Tmax"].asFloat()     ==  9999.0         );

        const Data extra = data["Extra"];
        const Data extrastrings = extra["SomeStrings"];

        CHECK( extra["SomeBoolean"].asBoolean()    == true  );
        CHECK( extra["AnotherBoolean"].asBoolean() == false );
        CHECK( extra["AnInteger"].asInteger()      == 123   );

        CHECK( extrastrings[0].asString() == "Hello" );
        CHECK( extrastrings[1].asString() == "Hallo" );

        CHECK_NOTHROW( data["Species"].with("Name", "Almandine") );
        CHECK_NOTHROW( data["Species"].with("Name", "Andradite") );

        CHECK( data["Species"].with("Name", "Almandine")["Name"].asString() == "Almandine" );
        CHECK( data["Species"].with("Name", "Andradite")["Name"].asString() == "Andradite" );

        CHECK_THROWS( data["Species"].with("Name", "Calcite") );
    }

    SECTION("Checking dumping of Data objects to YAML formatted strings")
    {
        const Data data = Data::parseYaml(yaml_testing_string);

        // CHECK( data.dumpYaml() == YAML::Dump(YAML::Load(yaml_testing_string)) );
        CHECK( data.dumpYaml() == "" );
    }

    SECTION("Checking dumping of Data objects to JSON formatted strings")
    {
        const Data data = Data::parseJson(json_testing_string);

        CHECK( data.dumpJson() == nlohmann::json::parse(json_testing_string).dump(2) ); // indent=2
    }

    SECTION("Checking encoding/decoding of custom types to/from Data objects")
    {
        const auto str = R"#(
            boolval: true
            doubleval: 1.2
            paramval: 3.4
            stringval: Moon
        )#";

        ComplicatedStuff stuff = Data::parse(str).as<ComplicatedStuff>();

        CHECK( stuff.boolval == true);
        CHECK( stuff.doubleval == 1.2);
        CHECK( stuff.paramval == 3.4);
        CHECK( stuff.stringval == "Moon");
        CHECK( stuff.vecval.empty() );
        CHECK( stuff.mapval.empty() );

        stuff.boolval = false;
        stuff.doubleval = 3.2;
        stuff.paramval = 7.6;
        stuff.stringval = "Sun";
        stuff.vecval = Strings{"A", "B"};
        stuff.mapval = Map<String, double>{{"U", 1.2}, {"V", 5.5}};

        Data encoded = stuff;

        CHECK( encoded.isDict() );
        CHECK( encoded.asDict().size() == 6 );

        CHECK( encoded.exists("boolval") );
        CHECK( encoded.exists("doubleval") );
        CHECK( encoded.exists("paramval") );
        CHECK( encoded.exists("stringval") );
        CHECK( encoded.exists("vecval") );
        CHECK( encoded.exists("mapval") );

        CHECK( encoded["boolval"].asBoolean()   == false );
        CHECK( encoded["doubleval"].asFloat()   == 3.2   );
        CHECK( encoded["paramval"].asParam()    == 7.6   );
        CHECK( encoded["stringval"].asString()  == "Sun" );
        CHECK( encoded["vecval"][0].asString()  == "A"   );
        CHECK( encoded["vecval"][1].asString()  == "B"   );
        CHECK( encoded["mapval"]["U"].asFloat() == 1.2   );
        CHECK( encoded["mapval"]["V"].asFloat() == 5.5   );
    }
}
