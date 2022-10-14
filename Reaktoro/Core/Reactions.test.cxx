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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Reactions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

namespace test {

class ReactionGeneratorUsingClass
{
public:
    auto operator()(SpeciesList const& species) const -> Vec<Reaction>
    {
        return { Reaction().withEquation("H2O(aq) = H+(aq) + OH-(aq)") };
    }
};

auto ReactionGeneratorUsingFunction(SpeciesList const& species) -> Vec<Reaction>
{
    return { Reaction().withEquation("H2O(aq) = H2(aq) + 0.5*O2(aq)") };
}

} // namespace test

TEST_CASE("Testing Reactions class", "[Reactions]")
{
    ChemicalSystem system = test::createChemicalSystem();
    Database db = system.database();

    Reactions reactions;

    reactions.add(db.reaction("NaCl(s) = Na+(aq) + Cl-(aq)"));
    reactions.add(db.reaction("CaCO3(s)"));
    reactions.add(test::ReactionGeneratorUsingClass());
    reactions.add(test::ReactionGeneratorUsingFunction);

    auto converted = reactions.convert(system.species());

    CHECK( converted.size() == 4 );
    CHECK( String(converted[0].equation()) == "NaCl(s) = Na+(aq) + Cl-(aq)" );
    CHECK( String(converted[1].equation()) == "CaCO3(s)" );
    CHECK( String(converted[2].equation()) == "H2O(aq) = H+(aq) + OH-(aq)" );
    CHECK( String(converted[3].equation()) == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );
}
