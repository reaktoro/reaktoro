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

TEST_CASE("Testing Core.yaml", "[Core.yaml]")
{
    std::string species_str1 = R"(
Name: Akermanite
Formula: Ca2MgSi2O7
Elements: {Ca: 2, Mg: 1, Si: 2, O: 7}
AggregateState: Solid
StandardThermoModel:
  MaierKelley:
    Gf: -3679250.6
    Hf: -3876463.4
    Sr: 209.32552
    Vr: 9.281e-05
    a: 251.41656
    b: 0.0476976
    c: -4769760.0
    Tmax: 1700.0
)";

    // Species species1 = yaml(species_str1);
}
