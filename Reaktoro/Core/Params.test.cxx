// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
    Params params;
    params.append(Param().id("A").value(1.0));
    params.append("B", 2.0);

    CHECK( params.size() == 2 );

    CHECK( params.find("A") == 0 );
    CHECK( params.find("B") == 1 );
    CHECK( params.find("C") == 2 );
    CHECK( params.find("D") == 2 );

    CHECK( params.index("A") == 0 );
    CHECK( params.index("B") == 1 );
    CHECK_THROWS( params.index("C") );
    CHECK_THROWS( params.index("D") );

    CHECK( params.get("A").value() == 1.0 );
    CHECK( params.get("B").value() == 2.0 );
    CHECK_THROWS( params.get("C") );
    CHECK_THROWS( params.get("D") );

    CHECK( params["A"].value() == 1.0 );
    CHECK( params["B"].value() == 2.0 );
    CHECK_THROWS( params["C"] );
    CHECK_THROWS( params["D"] );

    CHECK( params.exists("A") == true );
    CHECK( params.exists("B") == true );
    CHECK( params.exists("C") == false );
    CHECK( params.exists("D") == false );

    params = { Param("C", 4.0), 6.0 };

    CHECK( params.size() == 2 );
    CHECK( params[0].id() == "C" );
    CHECK( params[0].value() == 4.0 );
    CHECK( params[1].id() == "" );
    CHECK( params[1].value() == 6.0 );

    params.resize(3);
    CHECK( params.size() == 3 );
    CHECK( params[0].id() == "C" );
    CHECK( params[0].value() == 4.0 );
    CHECK( params[1].id() == "" );
    CHECK( params[1].value() == 6.0 );
    CHECK( params[2].id()    == ""  ); // the id of the newly added Param object!
    CHECK( params[2].value() == 0.0 ); // the value of the newly added Param object!

    CHECK( autodiff::detail::isVector<Params> );
}
