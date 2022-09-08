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
#include "StandardThermoModelFromData.hpp"
using namespace Reaktoro;

auto hkf = R"(
HKF:
  Gf: 1.0
  Hf: 2.0
  Sr: 0.0
  a1: 0.0
  a2: 0.0
  a3: 0.0
  a4: 0.0
  c1: 0.0
  c2: 0.0
  wref: 0.0
  charge: 0.0
  Tmax: 0.0
)";

auto mk = R"(
MaierKelley:
  Gf: 3.0
  Hf: 4.0
  Sr: 0.0
  Vr: 0.0
  a: 0.0
  b: 0.0
  c: 0.0
  Tmax: 0.0
)";

auto hp = R"(
HollandPowell:
  Gf: 5.0
  Hf: 6.0
  Sr: 0.0
  Vr: 0.0
  a: 0.0
  b: 0.0
  c: 0.0
  d: 0.0
  alpha0: 0.0
  kappa0: 0.0
  kappa0p: 0.0
  kappa0pp: 0.0
  numatoms: 0.0
  Tmax: 0.0
)";

auto mineral_hkf = R"(
MineralHKF:
  Gf: 7.0
  Hf: 8.0
  Sr: 0.0
  Vr: 0.0
  ntr: 3
  a: [0.0, 0.0, 0.0, 0.0]
  b: [0.0, 0.0, 0.0, 0.0]
  c: [0.0, 0.0, 0.0, 0.0]
  Ttr: [400, 500, 600]
  Htr: [0.0, 0.0, 0.0]
  Vtr: [0.0, 0.0, 0.0]
  dPdTtr: [0.0, 0.0, 0.0]
  Tmax: 0.0
)";

auto water_hkf = R"(
WaterHKF:
  Ttr: 273.16
  Str: 63.312288
  Gtr: -235517.36
  Htr: -287721.128
)";

auto non_existing_model = R"(
NonExistingModel:
  A1: 1.0
  A2: 2.0
  A3: 3.0
  A4: 4.0
)";

auto not_dict = R"(
- A
- B
)";

auto dict_but_not_single = R"(
A:
  A0: 1.0
B:
  B0: 2.0
)";

TEST_CASE("Testing StandardThermoModelFromData class", "[StandardThermoModelFromData]")
{
    const auto T = 25.0 + 273.15;
    const auto P =  1.0 * 1e5;

    WHEN("StandardThermoModel name is HKF")
    {
        auto model =  StandardThermoModelFromData(Data::parseYaml(hkf));
        auto props = model(T, P);

        CHECK(props.G0 == 1.0);
        CHECK(props.H0 == 2.0);
    }

    WHEN("StandardThermoModel name is MaierKelley")
    {
        auto model =  StandardThermoModelFromData(Data::parseYaml(mk));
        auto props = model(T, P);

        CHECK(props.G0 == 3.0);
        CHECK(props.H0 == 4.0);
    }

    WHEN("StandardThermoModel name is HollandPowell")
    {
        auto model =  StandardThermoModelFromData(Data::parseYaml(hp));
        auto props = model(T, P);

        CHECK(props.G0 == 5.0);
        CHECK(props.H0 == 6.0);
    }

    WHEN("StandardThermoModel name is MineralHKF")
    {
        auto model =  StandardThermoModelFromData(Data::parseYaml(mineral_hkf));
        auto props = model(T, P);

        CHECK(props.G0 == 7.0);
        CHECK(props.H0 == 8.0);
    }

    WHEN("StandardThermoModel name is WaterHKF")
    {
        auto model =  StandardThermoModelFromData(Data::parseYaml(water_hkf));
        auto props = model(T, P);

        CHECK(props.G0 == Approx(-237182));
        CHECK(props.H0 == Approx(-285831));
    }

    WHEN("yaml node is not valid")
    {
        CHECK_THROWS( StandardThermoModelFromData(Data::parseYaml(non_existing_model)) );
        CHECK_THROWS( StandardThermoModelFromData(Data::parseYaml(not_dict)) );
        CHECK_THROWS( StandardThermoModelFromData(Data::parseYaml(dict_but_not_single)) );
    }
}
