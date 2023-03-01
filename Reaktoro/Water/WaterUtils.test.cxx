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
#include <Reaktoro/Water/WaterUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing water utility methods", "[WaterUtils]")
{
    auto T = 273.15;
    auto P = 1e5;

    CHECK( waterDensityWagnerPruss(T, P, StateOfMatter::Gas) == Approx(0.0963975) );
    CHECK( waterDensityWagnerPruss(T + 25, P, StateOfMatter::Gas) == Approx(0.340528) );
    CHECK( waterDensityWagnerPruss(T + 50, P, StateOfMatter::Gas) == Approx(0.750385) );
    CHECK( waterDensityWagnerPruss(T + 75, P, StateOfMatter::Gas) == Approx(0.637419) );
    CHECK( waterDensityWagnerPruss(T + 100, P, StateOfMatter::Gas) == Approx(0.589669) );
    CHECK( waterDensityWagnerPruss(T + 400, P, StateOfMatter::Gas) == Approx(0.322301) );

    CHECK( waterDensityWagnerPruss(T, P, StateOfMatter::Liquid) == Approx(999.842) );
    CHECK( waterDensityWagnerPruss(T + 200, P, StateOfMatter::Liquid) == Approx(863.542) );
    CHECK( waterDensityWagnerPruss(T + 300, P, StateOfMatter::Liquid) == Approx(689.706) );
    CHECK( waterDensityWagnerPruss(T + 350, P, StateOfMatter::Liquid) == Approx(520.556) );
    CHECK( waterDensityWagnerPruss(T + 400, P, StateOfMatter::Liquid) == Approx(0.322301) );
    CHECK( waterDensityWagnerPruss(T + 500, P, StateOfMatter::Liquid) == Approx(0.280463) );
}
