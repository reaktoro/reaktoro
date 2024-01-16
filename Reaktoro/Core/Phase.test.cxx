// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Utils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Phase", "[Phase]")
{
    ActivityModel activity_model = [](ActivityPropsRef props, ActivityModelArgs args) {};

    Phase phase;
    phase = phase.withActivityModel(activity_model);

    REQUIRE( phase.activityModel() );

    SECTION("Testing PhasePhase::withSpecies with aqueous species")
    {
        phase = phase.withName("AqueousPhase");
        phase = phase.withSpecies(SpeciesList("H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--"));
        phase = phase.withStateOfMatter(StateOfMatter::Liquid);

        REQUIRE( phase.name() == "AqueousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.elements().size() == 3 );
        REQUIRE( phase.element(0).symbol() == "H" );
        REQUIRE( phase.element(1).symbol() == "C" );
        REQUIRE( phase.element(2).symbol() == "O" );
        REQUIRE( phase.species().size() == 8 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.species(3).name() == "H2(aq)"  );
        REQUIRE( phase.species(4).name() == "O2(aq)"  );
        REQUIRE( phase.species(5).name() == "CO2(aq)" );
        REQUIRE( phase.species(6).name() == "HCO3-"   );
        REQUIRE( phase.species(7).name() == "CO3--"   );
        REQUIRE( phase.speciesMolarMasses().isApprox(detail::molarMasses(phase.species())) );
    }

    SECTION("Testing PhasePhase::withSpecies with gaseous species")
    {
        Phase phase;
        phase = phase.withName("GaseousPhase");
        phase = phase.withSpecies(SpeciesList("H2O(g) CO2(g) CH4(g) O2(g) H2(g)"));
        phase = phase.withStateOfMatter(StateOfMatter::Gas);
        phase = phase.withActivityModel(activity_model);

        REQUIRE( phase.name() == "GaseousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Gas );
        REQUIRE( phase.elements().size() == 3 );
        REQUIRE( phase.element(0).symbol() == "H" );
        REQUIRE( phase.element(1).symbol() == "C" );
        REQUIRE( phase.element(2).symbol() == "O" );
        REQUIRE( phase.species().size() == 5 );
        REQUIRE( phase.species(0).name() == "H2O(g)" );
        REQUIRE( phase.species(1).name() == "CO2(g)" );
        REQUIRE( phase.species(2).name() == "CH4(g)" );
        REQUIRE( phase.species(3).name() == "O2(g)"  );
        REQUIRE( phase.species(4).name() == "H2(g)"  );
        REQUIRE( phase.speciesMolarMasses().isApprox(detail::molarMasses(phase.species())) );
    }

    SECTION("Testing PhasePhase::withSpecies with a single mineral species")
    {
        Phase phase;
        phase = phase.withName("Calcite");
        phase = phase.withSpecies({ Species("CaCO3(s)").withName("Calcite") });
        phase = phase.withStateOfMatter(StateOfMatter::Solid);

        REQUIRE( phase.name() == "Calcite" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phase.elements().size() == 3 );
        REQUIRE( phase.element(0).symbol() == "C" );
        REQUIRE( phase.element(1).symbol() == "O" );
        REQUIRE( phase.element(2).symbol() == "Ca" );
        REQUIRE( phase.species().size() == 1 );
        REQUIRE( phase.species(0).name() == "Calcite" );
        REQUIRE( phase.speciesMolarMasses().isApprox(detail::molarMasses(phase.species())) );
    }

    SECTION("Testing PhasePhase::withSpecies with species having uncommon aggregate states")
    {
        Phase phase;
        phase = phase.withName("AqueousPhase");
        phase = phase.withStateOfMatter(StateOfMatter::Liquid);

        REQUIRE_THROWS( phase.withSpecies(SpeciesList("H2O H+ OH-")) );
        // The aggregate state of H2O above cannot be inferred,
        // so it takes the value AggregateState::Undefined.
        // However, the aggregate state of the ionic species are
        // inferred as AggregateState::Aqueous by default.
        // Because of this mismatch in aggregate states, the above
        // Phase::withSpecies call should raise a runtime error.
    }
}
