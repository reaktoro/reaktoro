// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/Reaction.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Reaction class", "[Reaction]")
{
    Reaction reaction;

    reaction = reaction.withName("Dolomite");
    reaction = reaction.withEquation("CaCO3(s) = Ca++ + CO3--");
    reaction = reaction.withEquilibriumConstantFn([](real T, real P) { return T*P; });
    reaction = reaction.withRateFn([](const ChemicalProps& props) { return 1.0; });

    REQUIRE( reaction.name() == "Dolomite" );
    REQUIRE( reaction.equation().size() == 3 );
    REQUIRE( reaction.equation().coefficient("CaCO3(s)") == -1 );
    REQUIRE( reaction.equation().coefficient("Ca++") == 1 );
    REQUIRE( reaction.equation().coefficient("CO3--") == 1 );
    REQUIRE( reaction.equilibriumConstantFn() );
    REQUIRE( reaction.rateFn() );
}
