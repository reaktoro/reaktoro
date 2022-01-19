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
#include <Reaktoro/Core/Params.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Params class", "[Params]")
{
    Params foo;
    foo.set("A", 1.0);
    foo.set("B", 2.0);

    Params bar;
    bar.set("C", 5.0);
    bar.set("D", 6.0);
    bar.set("E", 7.0);

    Params params;
    params.set("Foo", foo);
    params.set("Bar", bar);

    CHECK( foo.size() == 2 );
    CHECK( foo.get("A").value() == 1.0 );
    CHECK( foo.get("B").value() == 2.0 );
    CHECK( foo.exists("A") == true );
    CHECK( foo.exists("B") == true );
    CHECK( foo.exists("C") == false );

    CHECK( bar.size() == 3 );
    CHECK( bar.get("C").value() == 5.0 );
    CHECK( bar.get("D").value() == 6.0 );
    CHECK( bar.get("E").value() == 7.0 );
    CHECK( bar.exists("C") == true );
    CHECK( bar.exists("D") == true );
    CHECK( bar.exists("E") == true );
    CHECK( bar.exists("F") == false );

    CHECK( params.size() == 5 );
    CHECK( params.at("Foo").get("A").value() == 1.0 );
    CHECK( params.at("Foo").get("B").value() == 2.0 );
    CHECK( params.at("Bar").get("C").value() == 5.0 );
    CHECK( params.at("Bar").get("D").value() == 6.0 );
    CHECK( params.at("Bar").get("E").value() == 7.0 );
    CHECK( params.exists("Foo") == true );
    CHECK( params.exists("Bar") == true );
    CHECK( params.exists("Joe") == false );
    CHECK( params.at("Foo").exists("A") == true );
    CHECK( params.at("Foo").exists("Z") == false );
    CHECK( params.at("Bar").exists("C") == true );
    CHECK( params.at("Bar").exists("Z") == false );
}
