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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ElementalComposition.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ElementalComposition", "[ElementalComposition]")
{
    ElementalComposition elements;

    //-----------------------------------------------------------------
    // Testing constructor ElementalComposition(Pairs<Element, double>)
    //-----------------------------------------------------------------
    elements = ElementalComposition({
        { Element("H"), 1.0 },
        { Element("C"), 2.0 },
        { Element("O"), 3.0 },
    });

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::symbols()
    //-----------------------------------------------------------------
    REQUIRE( identical(elements.symbols(), Strings{"H", "C", "O"}) );

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::coefficients()
    //-----------------------------------------------------------------
    REQUIRE( identical(elements.coefficients(), Vec<double>{1.0, 2.0, 3.0}) );

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::coefficient()
    //-----------------------------------------------------------------
    REQUIRE( elements.coefficient("H") == 1.0 );
    REQUIRE( elements.coefficient("C") == 2.0 );
    REQUIRE( elements.coefficient("O") == 3.0 );
    REQUIRE( elements.coefficient("X") == 0.0 );

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::repr()
    //-----------------------------------------------------------------
    REQUIRE( elements.repr() == "1:H 2:C 3:O" );

    //-----------------------------------------------------------------
    // Testing constructor ElementalComposition(Pairs<String, double>)
    //-----------------------------------------------------------------
    elements = ElementalComposition({
        { Element("Ca"), 4.0 },
        { Element("Mg"), 5.0 },
        { Element("F"),  6.0 },
    });

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::symbols()
    //-----------------------------------------------------------------
    REQUIRE( identical(elements.symbols(), Strings{"Ca", "Mg", "F"}) );

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::coefficients()
    //-----------------------------------------------------------------
    REQUIRE( identical(elements.coefficients(), Vec<double>{4.0, 5.0, 6.0}) );

    //-----------------------------------------------------------------
    // Testing method ElementalComposition::coefficient()
    //-----------------------------------------------------------------
    REQUIRE( elements.coefficient("Ca") == 4.0 );
    REQUIRE( elements.coefficient("Mg") == 5.0 );
    REQUIRE( elements.coefficient("F")  == 6.0 );
    REQUIRE( elements.coefficient("X")  == 0.0 );
}
