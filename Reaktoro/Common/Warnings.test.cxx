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
#include <Reaktoro/Common/Warnings.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Warnings", "[Warnings]")
{
    CHECK( Warnings::isEnabled(123) ); // by default, all warnings are enabled!

    Warnings::disable(123);
    CHECK( Warnings::isDisabled(123) );

    Warnings::enable(123);
    CHECK( Warnings::isEnabled(123) );

    Warnings::disable(124);
    Warnings::disable(179);

    // Check disabling already disabled warnings return false
    CHECK_FALSE( Warnings::disable(124) );
    CHECK_FALSE( Warnings::disable(179) );

    // Check enabling already enabled warnings return false
    CHECK_FALSE( Warnings::enable(897) );
    CHECK_FALSE( Warnings::enable(976) );
}
