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
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/Support/DatabaseParser.hpp>
using namespace Reaktoro;

namespace {

String doc = R"(
Elements:
  A:
    Symbol: A
    MolarMass: 1.0
  B:
    Symbol: B
    MolarMass: 2.0
Species:
  A2B(l):
    Name: A2B(l)
    Formula: A2B
    Elements: 2:A 1:B
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
  A2B3(aq):
    Name: A2B3(aq)
    Formula: A2B3
    Elements: 2:A 3:B
    AggregateState: Aqueous
    StandardThermoModel:
      HKF:
        Gf: 0.1
        Hf: 0.2
        Sr: 0.3
        a1: 0.4
        a2: 0.5
        a3: 0.6
        a4: 0.7
        c1: 0.8
        c2: 0.9
        wref: 1.0
        charge: 1.0
        Tmax: 1.2
  A2B(g):
    Name: A2B(g)
    Formula: A2B
    Elements: 2:A 1:B
    AggregateState: Gas
    FormationReaction:
      Reactants: 1:A2B(l)
      ReactionStandardThermoModel:
        ConstLgK:
          lgKr: 3.0
  A6B7(aq):
    Name: A6B7(aq)
    Formula: A6B7
    Elements: 6:A 7:B
    AggregateState: Aqueous
    FormationReaction:
      Reactants: 1:A2B(l) 2:A2B3(aq)
      ReactionStandardThermoModel:
        VantHoff:
          lgKr: 7.0
          dHr: 77.0
)";

} // namespace (anonymous)

TEST_CASE("Testing DatabaseParser class", "[DatabaseParser]")
{
    Data data = Data::fromYaml(doc);

    DatabaseParser db(data);

    auto elements = db.elements();
    auto species = db.species();

    CHECK( elements.size() == 2 );

    CHECK( elements[0].symbol() == "A" );
    CHECK( elements[0].molarMass() == 1.0 );
    CHECK( elements[1].symbol() == "B" );
    CHECK( elements[1].molarMass() == 2.0 );

    CHECK( species.size() == 4 );

    CHECK( species[0].name() == "A2B(l)" );
    CHECK( species[0].formula() == "A2B" );
    CHECK( species[0].elements().coefficient("A") == 2.0 );
    CHECK( species[0].elements().coefficient("B") == 1.0 );
    CHECK( species[0].aggregateState() == AggregateState::Liquid );
    CHECK( species[0].standardThermoModel().serialize()["MaierKelley"].isNull() );

    CHECK( species[1].name() == "A2B3(aq)" );
    CHECK( species[1].formula() == "A2B3" );
    CHECK( species[1].elements().coefficient("A") == 2.0 );
    CHECK( species[1].elements().coefficient("B") == 3.0 );
    CHECK( species[1].aggregateState() == AggregateState::Aqueous );
    CHECK( species[1].standardThermoModel().serialize()["HKF"].isNull() );

    CHECK( species[2].name() == "A2B(g)" );
    CHECK( species[2].formula() == "A2B" );
    CHECK( species[2].elements().coefficient("A") == 2.0 );
    CHECK( species[2].elements().coefficient("B") == 1.0 );
    CHECK( species[2].aggregateState() == AggregateState::Gas );
    CHECK( species[2].reaction().reactants().size() == 1 );
    CHECK( species[2].reaction().stoichiometry("A2B(l)") == 1 );

    CHECK( species[3].name() == "A6B7(aq)" );
    CHECK( species[3].formula() == "A6B7" );
    CHECK( species[3].elements().coefficient("A") == 6.0 );
    CHECK( species[3].elements().coefficient("B") == 7.0 );
    CHECK( species[3].aggregateState() == AggregateState::Aqueous );
    CHECK( species[3].reaction().reactants().size() == 2 );
    CHECK( species[3].reaction().stoichiometry("A2B(l)") == 1 );
    CHECK( species[3].reaction().stoichiometry("A2B3(aq)") == 2 );
}
