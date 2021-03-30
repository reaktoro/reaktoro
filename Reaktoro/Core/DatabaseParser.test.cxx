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

// Reaktoro includes
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/DatabaseParser.hpp>
using namespace Reaktoro;

std::string doc = R"(
Elements:
  - Symbol: A
    AtomicWeight: 1.0
  - Symbol: B
    AtomicWeight: 2.0
Species:
  - Name: A2B
    Formula: A2B
    Elements: {A: 2, B: 1}
    Charge: 0.0
    AggregateState: Liquid
    StandardThermoModel:
      MaierKelley:
        Gf: 1.234
        Hf: 2.345
        Sr: 0.0
        Vr: 0.0
        a: 0.0
        b: 0.0
        c: 0.0
        Tmax: 0.0
)";

TEST_CASE("Testing DatabaseParser class", "[DatabaseParser]")
{
    yaml node(doc);
    DatabaseParser dbparser(node);

}
