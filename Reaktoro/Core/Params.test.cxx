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
#include <Reaktoro/Core/Params.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Params class", "[Params]")
{
    Params foo;
    foo.append(Param().id("A").value(1.0));
    foo.append("B", 2.0);

    CHECK( foo.size() == 2 );

    CHECK( foo.find("A") == 0 );
    CHECK( foo.find("B") == 1 );
    CHECK( foo.find("C") == 2 );
    CHECK( foo.find("D") == 2 );

    CHECK( foo.index("A") == 0 );
    CHECK( foo.index("B") == 1 );
    CHECK_THROWS( foo.index("C") );
    CHECK_THROWS( foo.index("D") );

    CHECK( foo.get("A").value() == 1.0 );
    CHECK( foo.get("B").value() == 2.0 );
    CHECK_THROWS( foo.get("C") );
    CHECK_THROWS( foo.get("D") );

    CHECK( foo.exists("A") == true );
    CHECK( foo.exists("B") == true );
    CHECK( foo.exists("C") == false );
    CHECK( foo.exists("D") == false );
}
