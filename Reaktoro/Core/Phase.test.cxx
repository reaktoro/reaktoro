// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Phase", "[Phase]")
{
    StandardThermoPropsFn std_thermo_props_fn = [](real T, real P, const Species& species)
    {
        StandardThermoProps props;
        props.G0  = 1.0 * (T * P);
        props.H0  = 2.0 * (T * P);
        props.V0  = 3.0 * (T * P);
        props.Cp0 = 4.0 * (T * P);
        props.Cv0 = 5.0 * (T * P);
        return props;
    };

    ActivityPropsFn activity_props_fn = [](ActivityProps props, real T, real P, ArrayXrConstRef x)
    {
        props.Vex  = 1.0 * (T * P);
        props.VexT = 2.0 * (T * P);
        props.VexP = 3.0 * (T * P);
        props.Gex  = 4.0 * (T * P);
        props.Hex  = 5.0 * (T * P);
        props.Cpex = 6.0 * (T * P);
        props.Cvex = 7.0 * (T * P);
        props.ln_g = 8.0 * x;
        props.ln_a = 9.0 * x;
    };

    Phase phase;
    phase = phase.withStandardThermoPropsFn(std_thermo_props_fn);
    phase = phase.withActivityPropsFn(activity_props_fn);

    REQUIRE( phase.standardThermoPropsFn() );
    REQUIRE( phase.activityPropsFn() );

    SECTION("Testing PhasePhase::withSpecies with aqueous species")
    {
        phase = phase.withName("AqueousSolution");
        phase = phase.withSpecies(SpeciesList("H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--"));
        phase = phase.withStateOfMatter(StateOfMatter::Liquid);

        REQUIRE( phase.name() == "AqueousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 8 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.species(3).name() == "H2(aq)"  );
        REQUIRE( phase.species(4).name() == "O2(aq)"  );
        REQUIRE( phase.species(5).name() == "CO2(aq)" );
        REQUIRE( phase.species(6).name() == "HCO3-"   );
        REQUIRE( phase.species(7).name() == "CO3--"   );
    }

    SECTION("Testing PhasePhase::withSpecies with gaseous species")
    {
        Phase phase;
        phase = phase.withName("GaseousSolution");
        phase = phase.withSpecies(SpeciesList("H2O(g) CO2(g) CH4(g) O2(g) H2(g)"));
        phase = phase.withStateOfMatter(StateOfMatter::Gas);
        phase = phase.withStandardThermoPropsFn(std_thermo_props_fn);
        phase = phase.withActivityPropsFn(activity_props_fn);

        REQUIRE( phase.name() == "GaseousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Gas );
        REQUIRE( phase.species().size() == 5 );
        REQUIRE( phase.species(0).name() == "H2O(g)" );
        REQUIRE( phase.species(1).name() == "CO2(g)" );
        REQUIRE( phase.species(2).name() == "CH4(g)" );
        REQUIRE( phase.species(3).name() == "O2(g)"  );
        REQUIRE( phase.species(4).name() == "H2(g)"  );
    }

    SECTION("Testing PhasePhase::withSpecies with a single mineral species")
    {
        Phase phase;
        phase = phase.withName("Calcite");
        phase = phase.withSpecies({ Species("CaCO3(s)").withName("Calcite") });
        phase = phase.withStateOfMatter(StateOfMatter::Solid);

        REQUIRE( phase.name() == "Calcite" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phase.species().size() == 1 );
        REQUIRE( phase.species(0).name() == "Calcite" );
    }

    SECTION("Testing PhasePhase::withSpecies with species having uncommon aggregate states")
    {
        Phase phase;
        phase = phase.withName("AqueousSolution");
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
