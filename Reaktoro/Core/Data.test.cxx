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

const auto yaml_testing_string = R"(
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
)";

const auto json_testing_string = R"(
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
)";

TEST_CASE("Testing Data class", "[Data]")
{
    SECTION("testing construction and conversion methods")
    {
        Data foo;
        foo.add("A", 1.0);
        foo.add("B", 2);

        Data bar;
        bar.add("C", true);
        bar.add("D", "X");
        bar.add("E", Param(7.0));

        Data doo;
        doo.add(3.0);
        doo.add(6.0);
        doo.add(9.0);

        Data params;
        params.add("Foo", foo);
        params.add("Bar", bar);

        CHECK( foo["A"].number() == 1.0 );
        CHECK( foo["B"].integer() == 2 );
        CHECK( foo.exists("A") == true );
        CHECK( foo.exists("B") == true );
        CHECK( foo.exists("C") == false );

        CHECK( bar["C"].boolean() == true );
        CHECK( bar["D"].string() == "X" );
        CHECK( bar["E"].param()  == 7.0 );
        CHECK( bar.exists("C") == true );
        CHECK( bar.exists("D") == true );
        CHECK( bar.exists("E") == true );
        CHECK( bar.exists("F") == false );

        CHECK( doo[0].number() == 3.0 );
        CHECK( doo[1].number() == 6.0 );
        CHECK( doo[2].number() == 9.0 );

        CHECK( params["Foo"]["A"].number() == 1.0 );
        CHECK( params["Foo"]["B"].integer() == 2 );
        CHECK( params["Bar"]["C"].boolean() == true );
        CHECK( params["Bar"]["D"].string() == "X" );
        CHECK( params["Bar"]["E"].param()  == 7.0 );
        CHECK( params.exists("Foo") == true );
        CHECK( params.exists("Bar") == true );
        CHECK( params.exists("Joe") == false );
        CHECK( params["Foo"].exists("A") == true );
        CHECK( params["Foo"].exists("Z") == false );
        CHECK( params["Bar"].exists("C") == true );
        CHECK( params["Bar"].exists("Z") == false );
    }

    SECTION("testing construction using yaml and json strings")
    {
        const Data data = GENERATE(
            Data::fromYaml(yaml_testing_string),
            Data::fromJson(json_testing_string)
        );

        CHECK( data.isDict() );

        CHECK( data["Species"].isList() );
        CHECK( data["Species"].list().size() == 2 );

        const Data species0 = data["Species"][0];
        const Data species1 = data["Species"][1];

        CHECK( species0["Name"].string()           == "Almandine"           );
        CHECK( species0["Formula"].string()        == "Fe3Al2Si3O12"        );
        CHECK( species0["Elements"].string()       == "3:Fe 2:Al 3:Si 12:O" );
        CHECK( species0["AggregateState"].string() == "Solid"               );

        CHECK( species1["Name"].string()           == "Andradite"           );
        CHECK( species1["Formula"].string()        == "Ca3Fe2Si3O12"        );
        CHECK( species1["Elements"].string()       == "3:Ca 2:Fe 3:Si 12:O" );
        CHECK( species1["AggregateState"].string() == "Solid"               );

        const Data thermo0 = species0["StandardThermoModel"]["HollandPowell"];
        const Data thermo1 = species1["StandardThermoModel"]["HollandPowell"];

        CHECK( thermo0["Gf"].number()       == -4937500.0      );
        CHECK( thermo0["Hf"].number()       == -5260650.0      );
        CHECK( thermo0["Sr"].number()       ==  342.0          );
        CHECK( thermo0["Vr"].number()       ==  0.00011525     );
        CHECK( thermo0["a"].number()        ==  677.3          );
        CHECK( thermo0["b"].number()        ==  0.0            );
        CHECK( thermo0["c"].number()        == -3772700.0      );
        CHECK( thermo0["d"].number()        == -5044.0         );
        CHECK( thermo0["alpha0"].number()   ==  2.12e-05       );
        CHECK( thermo0["kappa0"].number()   ==  190000000000.0 );
        CHECK( thermo0["kappa0p"].number()  ==  2.98           );
        CHECK( thermo0["kappa0pp"].number() == -1.6e-11        );
        CHECK( thermo0["numatoms"].number() ==  20.0           );
        CHECK( thermo0["Tmax"].number()     ==  9999.0         );

        CHECK( thermo1["Gf"].number()       == -5426110.0      );
        CHECK( thermo1["Hf"].number()       == -5769080.0      );
        CHECK( thermo1["Sr"].number()       ==  316.4          );
        CHECK( thermo1["Vr"].number()       ==  0.00013204     );
        CHECK( thermo1["a"].number()        ==  638.6          );
        CHECK( thermo1["b"].number()        ==  0.0            );
        CHECK( thermo1["c"].number()        == -4955100.0      );
        CHECK( thermo1["d"].number()        == -3989.2         );
        CHECK( thermo1["alpha0"].number()   ==  2.86e-05       );
        CHECK( thermo1["kappa0"].number()   ==  158800000000.0 );
        CHECK( thermo1["kappa0p"].number()  ==  5.68           );
        CHECK( thermo1["kappa0pp"].number() == -3.6e-11        );
        CHECK( thermo1["numatoms"].number() ==  20.0           );
        CHECK( thermo1["Tmax"].number()     ==  9999.0         );

        const Data extra = data["Extra"];
        const Data extrastrings = extra["SomeStrings"];

        CHECK( extra["SomeBoolean"].boolean()    == true  );
        CHECK( extra["AnotherBoolean"].boolean() == false );
        CHECK( extra["AnInteger"].integer()      == 123   );

        CHECK( extrastrings[0].string() == "Hello" );
        CHECK( extrastrings[1].string() == "Hallo" );

        CHECK_NOTHROW( data["Species"].with("Name", "Almandine") );
        CHECK_NOTHROW( data["Species"].with("Name", "Andradite") );

        CHECK( data["Species"].with("Name", "Almandine")["Name"].string() == "Almandine" );
        CHECK( data["Species"].with("Name", "Andradite")["Name"].string() == "Andradite" );

        CHECK_THROWS( data["Species"].with("Name", "Calcite") );
    }
}
