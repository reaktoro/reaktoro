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
#include <Reaktoro/Core/Param.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Param class", "[Param]")
{
    Param x;

    CHECK( x.value() == 0.0 );
    CHECK( x.id() == "" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == false );

    x = Param("x", 1.0);

    CHECK( x.value() == 1.0 );
    CHECK( x.id() == "x" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == false );

    x = Param().value(3.0).id("xx").lowerbound(1.0).upperbound(7.0).isconst(true);

    CHECK( x.value() == 3.0 );
    CHECK( x.id() == "xx" );
    CHECK( x.lowerbound() == 1.0 );
    CHECK( x.upperbound() == 7.0 );
    CHECK( x.isconst() == true );

    x.lowerbound(0.0).id("y").upperbound(17.0).value(9.0).isconst(false);

    CHECK( x.value() == 9.0 );
    CHECK( x.id() == "y" );
    CHECK( x.lowerbound() == 0.0 );
    CHECK( x.upperbound() == 17.0 );
    CHECK( x.isconst() == false );

    Param xclone = x.clone();

    CHECK( xclone.value() == x.value() );
    CHECK( xclone.id() == x.id() );
    CHECK( xclone.lowerbound() == x.lowerbound() );
    CHECK( xclone.upperbound() == x.upperbound() );
    CHECK( xclone.isconst() == x.isconst() );

    x.value(16.0);
    xclone.value(13.0);

    CHECK( x.value() == 16.0 ); // changing xclone does not change x since they have independent underlying data!
    CHECK( xclone.value() == 13.0 );

    Param z = x; // z points to x

    CHECK( z.value() == 16.0 );
    CHECK( z.id() == "y" );
    CHECK( z.lowerbound() == 0.0 );
    CHECK( z.upperbound() == 17.0 );
    CHECK( z.isconst() == false );

    z.value(10.0).id("zzz"); // changing Param z here actually changes underlying Param x

    CHECK( x.value() == 10.0 ); // x.value() changed with z.value(10.0)
    CHECK( x.id() == "zzz" );   // x.id() changed with z.id("z")

    x = Param::Constant(11.0);

    CHECK( x.value() == 11.0 );
    CHECK( x.id() == "" );
    CHECK( x.lowerbound() < 0.0 );
    CHECK( x.upperbound() > 0.0 );
    CHECK( std::isinf(x.lowerbound()) );
    CHECK( std::isinf(x.upperbound()) );
    CHECK( x.isconst() == true );

    CHECK( autodiff::detail::Order<Param> == 1 );
}
