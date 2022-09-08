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

const auto some_model_params_1 = R"(
SomeModelParams1:
  A: 1.1
  B: 1.2
)";

const auto some_model_params_2 = R"(
SomeModelParams2:
  C: 2.1
  D: 2.2
)";

TEST_CASE("Testing Params class", "[Params]")
{
    Params params1(Data::fromYaml(some_model_params_1));
    Params params2(Data::fromYaml(some_model_params_2));

    Params params;
    params += params1;
    params += params2;

    CHECK( params["SomeModelParams1"]["A"].number() == 1.1 );
    CHECK( params["SomeModelParams1"]["B"].number() == 1.2 );

    CHECK( params["SomeModelParams2"]["C"].number() == 2.1 );
    CHECK( params["SomeModelParams2"]["D"].number() == 2.2 );
}
