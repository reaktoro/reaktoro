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
#include <Reaktoro/Core/Data.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Data class", "[Data]")
{
    Data foo;
    foo.add("A", 1.0);
    foo.add("B", 2);

    Data bar;
    bar.add("C", true);
    bar.add("D", "X");
    bar.add("E", Param(7.0));

    Data doo;
    doo.add(3.0);
    doo.add(6.0);
    doo.add(9.0);

    Data params;
    params.add("Foo", foo);
    params.add("Bar", bar);

    CHECK( foo["A"].number() == 1.0 );
    CHECK( foo["B"].integer() == 2 );
    CHECK( foo.exists("A") == true );
    CHECK( foo.exists("B") == true );
    CHECK( foo.exists("C") == false );

    CHECK( bar["C"].boolean() == true );
    CHECK( bar["D"].string() == "X" );
    CHECK( bar["E"].param()  == 7.0 );
    CHECK( bar.exists("C") == true );
    CHECK( bar.exists("D") == true );
    CHECK( bar.exists("E") == true );
    CHECK( bar.exists("F") == false );

    CHECK( doo[0].number() == 3.0 );
    CHECK( doo[1].number() == 6.0 );
    CHECK( doo[2].number() == 9.0 );

    CHECK( params["Foo"]["A"].number() == 1.0 );
    CHECK( params["Foo"]["B"].integer() == 2 );
    CHECK( params["Bar"]["C"].boolean() == true );
    CHECK( params["Bar"]["D"].string() == "X" );
    CHECK( params["Bar"]["E"].param()  == 7.0 );
    CHECK( params.exists("Foo") == true );
    CHECK( params.exists("Bar") == true );
    CHECK( params.exists("Joe") == false );
    CHECK( params["Foo"].exists("A") == true );
    CHECK( params["Foo"].exists("Z") == false );
    CHECK( params["Bar"].exists("C") == true );
    CHECK( params["Bar"].exists("Z") == false );
}
