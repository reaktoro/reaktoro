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
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
using namespace Reaktoro;

TEST_CASE("Testing EquilibriumResult", "[EquilibriumResult]")
{
    EquilibriumResult result;

    CHECK( result.succeeded() == false );
    CHECK( result.failed() == true );
    CHECK( result.iterations() == 0 );

    result.optima.succeeded = true;
    result.optima.iterations = 23;

    CHECK( result.succeeded() == true );
    CHECK( result.failed() == false );
    CHECK( result.iterations() == 23 );
}
