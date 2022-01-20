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
#include <Reaktoro/Core/ReactionEquation.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ReactionEquation class", "[ReactionEquation]")
{
    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ReactionEquation()
    //-------------------------------------------------------------------------
    ReactionEquation equation;

    REQUIRE( equation.empty() );
    REQUIRE( equation.size() == 0 );

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ReactionEquation(Pairs<Species, double>)
    //-------------------------------------------------------------------------
    equation = ReactionEquation({ {Species("H2O"), -1}, {Species("H+"), 1}, {Species("OH-"), 1}});

    REQUIRE( equation.size() == 3 );
    REQUIRE( equation.coefficient("H2O")  == -1.0 );
    REQUIRE( equation.coefficient("H+")   ==  1.0 );
    REQUIRE( equation.coefficient("OH-")  ==  1.0 );
    REQUIRE( equation.coefficient("Fe++") ==  0.0 );

    equation = ReactionEquation({ {Species("CaCl2"), -1}, {Species("Ca++"), 1}, {Species("Cl-"), 2}});

    REQUIRE( equation.size() == 3 );
    REQUIRE( equation.coefficient("CaCl2") == -1.0 );
    REQUIRE( equation.coefficient("Ca++")  ==  1.0 );
    REQUIRE( equation.coefficient("Cl-")   ==  2.0 );
    REQUIRE( equation.coefficient("Fe++")  ==  0.0 );

    equation = ReactionEquation({ {Species("Ca++"), -1}, {Species("Mg++"), -1}, {Species("HCO3-"), -2}, {Species("CaMg(CO3)2"), 1}, {Species("H+"), 2}});

    REQUIRE( equation.size() == 5 );
    REQUIRE( equation.coefficient("Ca++")       == -1.0 );
    REQUIRE( equation.coefficient("Mg++")       == -1.0 );
    REQUIRE( equation.coefficient("HCO3-")      == -2.0 );
    REQUIRE( equation.coefficient("CaMg(CO3)2") ==  1.0 );
    REQUIRE( equation.coefficient("H+")         ==  2.0 );
    REQUIRE( equation.coefficient("Fe++")       ==  0.0 );

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ReactionEquation(String)
    //-------------------------------------------------------------------------
    equation = ReactionEquation("H2O = H+ + OH-");

    REQUIRE( equation.size() == 3 );
    REQUIRE( equation.coefficient("H2O")  == -1.0 );
    REQUIRE( equation.coefficient("H+")   ==  1.0 );
    REQUIRE( equation.coefficient("OH-")  ==  1.0 );
    REQUIRE( equation.coefficient("Fe++") ==  0.0 );

    equation = ReactionEquation("CaCl2 = Ca++ + 2*Cl-");

    REQUIRE( equation.size() == 3 );
    REQUIRE( equation.coefficient("CaCl2") == -1.0 );
    REQUIRE( equation.coefficient("Ca++")  ==  1.0 );
    REQUIRE( equation.coefficient("Cl-")   ==  2.0 );
    REQUIRE( equation.coefficient("Fe++")  ==  0.0 );

    equation = ReactionEquation("Ca++ + Mg++ + 2*HCO3- = CaMg(CO3)2 + 2*H+");

    REQUIRE( equation.size() == 5 );
    REQUIRE( equation.coefficient("Ca++")       == -1.0 );
    REQUIRE( equation.coefficient("Mg++")       == -1.0 );
    REQUIRE( equation.coefficient("HCO3-")      == -2.0 );
    REQUIRE( equation.coefficient("CaMg(CO3)2") ==  1.0 );
    REQUIRE( equation.coefficient("Fe++")       ==  0.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: operator==(ReactionEquation,ReactionEquation)
    //-------------------------------------------------------------------------
    REQUIRE( ReactionEquation("H2O = H+ + OH-") == ReactionEquation("1*H2O = 1*H+ + 1*OH-") );
    REQUIRE( ReactionEquation("Ca++ + Mg++ + 2*HCO3- = CaMg(CO3)2 + 2*H+") == ReactionEquation("2*HCO3- + Mg++ + Ca++ = 2*H+ + CaMg(CO3)2") );

    REQUIRE_FALSE( ReactionEquation("Ca++ + Mg++ + 2*HCO3- = CaMg(CO3)2 + 2*H+") == ReactionEquation("CaMg(CO3)2 + 2*H+ = Ca++ + Mg++ + 2*HCO3-") );
    REQUIRE_FALSE( ReactionEquation("H2O = H+ + OH-") == ReactionEquation("H+ + OH- = H2O") );
}
