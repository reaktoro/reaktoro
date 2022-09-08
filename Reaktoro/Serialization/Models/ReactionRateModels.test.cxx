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
using namespace Catch;

// Reaktoro includes
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Serialization/Models/ReactionRateModels.hpp>
using namespace Reaktoro;

TEST_CASE("Testing serialization of ReactionRateModelParamsPalandriKharaka", "[Serialization][Models][ReactionRateModels]")
{
    ReactionRateModelParamsPalandriKharaka params;

    WHEN("there is only one mechanism")
    {
        String yml = R"#(
            Mineral: Halite
            Mechanisms:
              Neutral: { lgk: -0.21, E: 7.4 }
        )#";

        params = Data::parse(yml).as<ReactionRateModelParamsPalandriKharaka>();

        CHECK( params.mineral == "Halite" );
        CHECK( params.othernames.empty() );
        CHECK( params.mechanisms.size() == 1 );
        CHECK( params.mechanisms[0].name == "Neutral" );
        CHECK( params.mechanisms[0].lgk == -0.21 );
        CHECK( params.mechanisms[0].E == 7.40 );
        CHECK( params.mechanisms[0].p == 1.00 );
        CHECK( params.mechanisms[0].q == 1.00 );
        CHECK( params.mechanisms[0].catalysts.size() == 0 );

    }

    WHEN("there are two mechanisms")
    {
        String yml = R"#(
            Mineral: Pyrite
            Mechanisms:
              Acid:    { lgk: -7.52, E: 56.9, a(H+): -0.500, a(Fe+3): 0.500 }
              Neutral: { lgk: -4.55, E: 56.9, a(O2): 0.500 }
        )#";

        params = Data::parse(yml).as<ReactionRateModelParamsPalandriKharaka>();

        CHECK( params.mineral == "Pyrite" );
        CHECK( params.othernames.empty() );
        CHECK( params.mechanisms.size() == 2 );
        CHECK( params.mechanisms[0].name == "Acid" );
        CHECK( params.mechanisms[0].lgk == -7.52 );
        CHECK( params.mechanisms[0].E == 56.9 );
        CHECK( params.mechanisms[0].p == 1.00 );
        CHECK( params.mechanisms[0].q == 1.00 );
        CHECK( params.mechanisms[0].catalysts.size() == 2 );
        CHECK( params.mechanisms[0].catalysts[0].formula == "H+" );
        CHECK( params.mechanisms[0].catalysts[0].power == -0.50 );
        CHECK( params.mechanisms[0].catalysts[0].property == "a" );
        CHECK( params.mechanisms[0].catalysts[1].formula == "Fe+3" );
        CHECK( params.mechanisms[0].catalysts[1].power == 0.50 );
        CHECK( params.mechanisms[0].catalysts[1].property == "a" );
        CHECK( params.mechanisms[1].name == "Neutral" );
        CHECK( params.mechanisms[1].lgk == -4.55 );
        CHECK( params.mechanisms[1].E == 56.9 );
        CHECK( params.mechanisms[1].p == 1.00 );
        CHECK( params.mechanisms[1].q == 1.00 );
        CHECK( params.mechanisms[1].catalysts.size() == 1 );
        CHECK( params.mechanisms[1].catalysts[0].formula == "O2" );
        CHECK( params.mechanisms[1].catalysts[0].power == 0.50 );
        CHECK( params.mechanisms[1].catalysts[0].property == "a" );
    }

    WHEN("there are three mechanisms")
    {
        String yml = R"#(
            Mineral: Epidote
            Mechanisms:
              Acid:    { lgk: -10.60, E: 71.1, a(H+):  0.338 }
              Neutral: { lgk: -11.99, E: 70.7 }
              Base:    { lgk: -17.33, E: 79.1, a(H+): -0.556 }
        )#";

        params = Data::parse(yml).as<ReactionRateModelParamsPalandriKharaka>();

        CHECK( params.mineral == "Epidote" );
        CHECK( params.othernames.empty() );
        CHECK( params.mechanisms.size() == 3 );
        CHECK( params.mechanisms[0].name == "Acid" );
        CHECK( params.mechanisms[0].lgk == -10.60 );
        CHECK( params.mechanisms[0].E == 71.1 );
        CHECK( params.mechanisms[0].p == 1.00 );
        CHECK( params.mechanisms[0].q == 1.00 );
        CHECK( params.mechanisms[0].catalysts.size() == 1 );
        CHECK( params.mechanisms[0].catalysts[0].formula == "H+" );
        CHECK( params.mechanisms[0].catalysts[0].power == 0.338 );
        CHECK( params.mechanisms[0].catalysts[0].property == "a" );
        CHECK( params.mechanisms[1].name == "Neutral" );
        CHECK( params.mechanisms[1].lgk == -11.99 );
        CHECK( params.mechanisms[1].E == 70.7 );
        CHECK( params.mechanisms[1].p == 1.00 );
        CHECK( params.mechanisms[1].q == 1.00 );
        CHECK( params.mechanisms[1].catalysts.size() == 0 );
        CHECK( params.mechanisms[2].name == "Base" );
        CHECK( params.mechanisms[2].lgk == -17.33 );
        CHECK( params.mechanisms[2].E == 79.1 );
        CHECK( params.mechanisms[2].p == 1.00 );
        CHECK( params.mechanisms[2].q == 1.00 );
        CHECK( params.mechanisms[2].catalysts.size() == 1 );
        CHECK( params.mechanisms[2].catalysts[0].formula == "H+" );
        CHECK( params.mechanisms[2].catalysts[0].power == -0.556 );
        CHECK( params.mechanisms[2].catalysts[0].property == "a" );
    }

    WHEN("there are many mechanisms")
    {
        String yml = R"#(
            Mineral: FakeMineral
            OtherNames:
              - FakeMineral1
              - FakeMineral2
            Mechanisms:
              Acid:      { lgk: -0.30, E: 14.4, a(H+): 1.000 }
              Neutral:   { lgk: -5.81, E: 23.5 }
              Carbonate: { lgk: -3.48, E: 35.4, P(CO2): 1.000 }
              Base:      { lgk: -17.33, E: 79.1, a((Ca)(CO3)): -0.556, p: 2.34, q: 1.46 }
        )#";

        params = Data::parse(yml).as<ReactionRateModelParamsPalandriKharaka>();

        CHECK( params.mineral == "FakeMineral" );
        CHECK( params.othernames == Strings{"FakeMineral1", "FakeMineral2"} );
        CHECK( params.mechanisms.size() == 4 );
        CHECK( params.mechanisms[0].name == "Acid" );
        CHECK( params.mechanisms[0].lgk == -0.30 );
        CHECK( params.mechanisms[0].E == 14.4 );
        CHECK( params.mechanisms[0].p == 1.00 );
        CHECK( params.mechanisms[0].q == 1.00 );
        CHECK( params.mechanisms[0].catalysts.size() == 1 );
        CHECK( params.mechanisms[0].catalysts[0].formula == "H+" );
        CHECK( params.mechanisms[0].catalysts[0].power == 1.0 );
        CHECK( params.mechanisms[0].catalysts[0].property == "a" );
        CHECK( params.mechanisms[1].name == "Neutral" );
        CHECK( params.mechanisms[1].lgk == -5.81 );
        CHECK( params.mechanisms[1].E == 23.5 );
        CHECK( params.mechanisms[1].p == 1.00 );
        CHECK( params.mechanisms[1].q == 1.00 );
        CHECK( params.mechanisms[1].catalysts.size() == 0 );
        CHECK( params.mechanisms[2].name == "Carbonate" );
        CHECK( params.mechanisms[2].lgk == -3.48 );
        CHECK( params.mechanisms[2].E == 35.4 );
        CHECK( params.mechanisms[2].p == 1.00 );
        CHECK( params.mechanisms[2].q == 1.00 );
        CHECK( params.mechanisms[2].catalysts.size() == 1 );
        CHECK( params.mechanisms[2].catalysts[0].formula == "CO2" );
        CHECK( params.mechanisms[2].catalysts[0].power == 1.0 );
        CHECK( params.mechanisms[2].catalysts[0].property == "P" );
        CHECK( params.mechanisms[3].name == "Base" );
        CHECK( params.mechanisms[3].lgk == -17.33 );
        CHECK( params.mechanisms[3].E == 79.1 );
        CHECK( params.mechanisms[3].p == 2.34 );
        CHECK( params.mechanisms[3].q == 1.46 );
        CHECK( params.mechanisms[3].catalysts.size() == 1 );
        CHECK( params.mechanisms[3].catalysts[0].formula == "(Ca)(CO3)" );
        CHECK( params.mechanisms[3].catalysts[0].power == -0.556 );
        CHECK( params.mechanisms[3].catalysts[0].property == "a" );
    }
}
