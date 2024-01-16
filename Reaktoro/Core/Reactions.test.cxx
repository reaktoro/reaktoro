// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Reactions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

namespace test {

auto createReactionRateModel(real value) -> ReactionRateModel
{
    return [=](ChemicalProps const& props) { return value; };
}

auto createReactionRateModelGenerator(real value) -> ReactionRateModelGenerator
{
    return [=](ReactionRateModelGeneratorArgs args) {
        return [=](ChemicalProps const& props) { return value; };
    };
}

class ReactionGeneratorUsingClass
{
public:
    auto operator()(ReactionGeneratorArgs args) const -> Vec<Reaction>
    {
        return { Reaction()
            .withEquation("H2O(aq) = H+(aq) + OH-(aq)")
            .withRateModel(createReactionRateModel(5.0))
        };
    }
};

auto ReactionGeneratorUsingFunction(ReactionGeneratorArgs args) -> Vec<Reaction>
{
    return { Reaction()
        .withEquation("H2O(aq) = H2(aq) + 0.5*O2(aq)")
        .withRateModel(createReactionRateModel(6.0))
    };
}

} // namespace test

TEST_CASE("Testing Reactions class", "[Reactions]")
{
    ChemicalSystem system = test::createChemicalSystem();
    ChemicalProps props(system);

    Database db = system.database();

    Reaction reaction1 = db.reaction("NaCl(s) = Na+(aq) + Cl-(aq)");
    Reaction reaction2 = db.reaction("CaCO3(s)");

    GeneralReaction generalreaction1("CO2(g) = CO2(aq)");
    GeneralReaction generalreaction2("HCO3-(aq) + H+(aq) = CO2(aq) + H2O(aq)");

    reaction1 = reaction1.withRateModel(test::createReactionRateModel(1.0));
    reaction2 = reaction2.withRateModel(test::createReactionRateModel(2.0));
    generalreaction1.setRateModel(test::createReactionRateModel(3.0));
    generalreaction2.setRateModel(test::createReactionRateModelGenerator(4.0));

    Reactions reactionsA;
    reactionsA.add(reaction1);
    reactionsA.add(reaction2);
    reactionsA.add(generalreaction1);
    reactionsA.add(generalreaction2);
    reactionsA.add(test::ReactionGeneratorUsingClass());
    reactionsA.add(test::ReactionGeneratorUsingFunction);

    Reactions reactionsB(
        reaction1,
        reaction2,
        generalreaction1,
        generalreaction2,
        test::ReactionGeneratorUsingClass(),
        test::ReactionGeneratorUsingFunction
    );

    auto checkReactionsConversion = [&](Reactions const& reactions)
    {
        ReactionGeneratorArgs args{ system.database(), system.species(), system.phases(), system.surfaces() };

        auto converted = reactions.convert(args);

        CHECK( converted.size() == 6 );

        CHECK( String(converted[0].equation()) == "NaCl(s) = Na+(aq) + Cl-(aq)" );
        CHECK( String(converted[1].equation()) == "CaCO3(s)" );
        CHECK( String(converted[2].equation()) == "CO2(g) = CO2(aq)" );
        CHECK( String(converted[3].equation()) == "HCO3-(aq) + H+(aq) = CO2(aq) + H2O(aq)" );
        CHECK( String(converted[4].equation()) == "H2O(aq) = H+(aq) + OH-(aq)" );
        CHECK( String(converted[5].equation()) == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );

        CHECK( converted[0].rate(props).val() == 1.0 );
        CHECK( converted[1].rate(props).val() == 2.0 );
        CHECK( converted[2].rate(props).val() == 3.0 );
        CHECK( converted[3].rate(props).val() == 4.0 );
        CHECK( converted[4].rate(props).val() == 5.0 );
        CHECK( converted[5].rate(props).val() == 6.0 );
    };

    checkReactionsConversion(reactionsA);
    checkReactionsConversion(reactionsB);
}
