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
#include <Reaktoro/Core/Embedded.hpp>
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

const auto some_overlapping_and_new_model_params = R"(
SomeModelParams1:
  B: 7.7
SomeModelParams2:
  C: 9.9
SomeModelParams3:
  E: 100.0
  F: 200.0
  G: [111.0, 222.0, 333.0]
)";

TEST_CASE("Testing Params class", "[Params]")
{
    SECTION("Check basic usage of a Params object")
    {
        Params params1(Data::parse(some_model_params_1));
        Params params2(Data::parse(some_model_params_2));
        Params params3(Data::parse(some_overlapping_and_new_model_params));

        Params params;
        params += params1;
        params += params2;

        CHECK( params["SomeModelParams1"]["A"].asFloat() == 1.1 );
        CHECK( params["SomeModelParams1"]["B"].asFloat() == 1.2 );
        CHECK( params["SomeModelParams2"]["C"].asFloat() == 2.1 );
        CHECK( params["SomeModelParams2"]["D"].asFloat() == 2.2 );

        params += params3;

        CHECK( params["SomeModelParams1"]["A"].asFloat()    == 1.1 );   // same value as before
        CHECK( params["SomeModelParams1"]["B"].asFloat()    == 7.7 );   // overwritten value!
        CHECK( params["SomeModelParams2"]["C"].asFloat()    == 9.9 );   // overwritten value!
        CHECK( params["SomeModelParams2"]["D"].asFloat()    == 2.2 );   // same value as before
        CHECK( params["SomeModelParams3"]["E"].asFloat()    == 100.0 ); // new inserted parameter!
        CHECK( params["SomeModelParams3"]["F"].asFloat()    == 200.0 ); // new inserted parameter!
        CHECK( params["SomeModelParams3"]["G"][0].asFloat() == 111.0 ); // new inserted parameter!
        CHECK( params["SomeModelParams3"]["G"][1].asFloat() == 222.0 ); // new inserted parameter!
        CHECK( params["SomeModelParams3"]["G"][2].asFloat() == 333.0 ); // new inserted parameter!
    }

    SECTION("Check when constructed from an embedded or local resource")
    {
        Params params = GENERATE(
            Params::embedded("PalandriKharaka.yaml"),
            Params::embedded("PalandriKharaka.json"),
            Params::local(REAKTORO_PARAMS_DIR"/PalandriKharaka.yaml"),
            Params::local(REAKTORO_PARAMS_DIR"/PalandriKharaka.json")
        );

        auto const& data = params.data();

        CHECK( data.exists("ReactionRateModelParams") );
        CHECK( data.at("ReactionRateModelParams").exists("PalandriKharaka") );
        CHECK( data.at("ReactionRateModelParams").at("PalandriKharaka").exists("Albite") );
    }
}
