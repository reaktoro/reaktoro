// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Common/StringUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StringUtils", "[StringUtils]")
{
    REQUIRE( stringfy(".", "a", "b") == "a.b");
    REQUIRE( stringfy(".", 1, 2, "abc") == "1.2.abc");
    REQUIRE( str(1, 2, " abc") == "12 abc");
}
