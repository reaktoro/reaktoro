// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
    // Testing constructor ElementalComposition(Map<Element, double>)
    //-----------------------------------------------------------------
    elements = ElementalComposition({
        { Element("H"), 1.0 },
        { Element("C"), 2.0 },
        { Element("O"), 3.0 },
    });

    // Testing method ElementalComposition::symbols()
    REQUIRE( elements.symbols().size() == 3 );
    REQUIRE( elements.symbols()[0] == "H" );
    REQUIRE( elements.symbols()[1] == "C" );
    REQUIRE( elements.symbols()[2] == "O" );

    // Testing method ElementalComposition::coefficients()
    REQUIRE( elements.coefficients().size() == 3 );
    REQUIRE( elements.coefficients()[0] == 1.0 );
    REQUIRE( elements.coefficients()[1] == 2.0 );
    REQUIRE( elements.coefficients()[2] == 3.0 );

    REQUIRE( contains(elements.symbols(), "C") );
    REQUIRE( contains(elements.symbols(), "O") );

    // Testing method ElementalComposition::coefficient()
    REQUIRE( elements.coefficient("H") == 1.0 );
    REQUIRE( elements.coefficient("C") == 2.0 );
    REQUIRE( elements.coefficient("O") == 3.0 );

    //-----------------------------------------------------------------
    // Testing constructor ElementalComposition(Map<String, double>)
    //-----------------------------------------------------------------
    elements = ElementalComposition(Map<String, double>{
        { "H", 1.0 },
        { "C", 2.0 },
        { "O", 3.0 },
    });

    // Testing method ElementalComposition::symbols()
    REQUIRE( elements.symbols().size() == 3 );
    REQUIRE( elements.symbols()[0] == "H" );
    REQUIRE( elements.symbols()[1] == "C" );
    REQUIRE( elements.symbols()[2] == "O" );

    // Testing method ElementalComposition::coefficients()
    REQUIRE( elements.coefficients().size() == 3 );
    REQUIRE( elements.coefficients()[0] == 1.0 );
    REQUIRE( elements.coefficients()[1] == 2.0 );
    REQUIRE( elements.coefficients()[2] == 3.0 );

    REQUIRE( contains(elements.symbols(), "C") );
    REQUIRE( contains(elements.symbols(), "O") );

    // Testing method ElementalComposition::coefficient()
    REQUIRE( elements.coefficient("H") == 1.0 );
    REQUIRE( elements.coefficient("C") == 2.0 );
    REQUIRE( elements.coefficient("O") == 3.0 );
}
