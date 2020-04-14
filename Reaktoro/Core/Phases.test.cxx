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
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Phases", "[Phases]")
{
    StandardThermoPropsFn std_thermo_props_fn = [](real T, real P, const Species& species)
    {
        return StandardThermoProps{ 1.0, 2.0, 3.0, 4.0, 5.0 }; // G0, H0, V0, Cp0, Cv0
    };

    ActivityModel activitymodel = [](const SpeciesList& species)
    {
        ActivityPropsFn fn = [](ActivityProps props, real T, real P, ArrayXrConstRef x)
        {
            props.Vex  = 123.0;
            props.VexT = 123.0;
            props.VexP = 123.0;
            props.Gex  = 123.0;
            props.Hex  = 123.0;
            props.Cpex = 123.0;
            props.Cvex = 123.0;
            props.ln_g = 123.0;
            props.ln_a = 123.0;
        };

        return fn;
    };

    Database db;
    db.addSpecies(Species("H2O(aq)"));
    db.addSpecies(Species("H+"));
    db.addSpecies(Species("OH-"));
    db.addSpecies(Species("Na+"));
    db.addSpecies(Species("Cl-"));

    ThermoEngine engine(db, std_thermo_props_fn);

    GenericPhase genericphase("H2O(aq) H+ OH-");
    genericphase.named("AqueousPhase");
    genericphase.set(StateOfMatter::Liquid);
    genericphase.set(AggregateState::Aqueous);
    genericphase.set(activitymodel);

    Phase phase = genericphase.convert(engine, {"H", "O", "Na", "Cl"});

    REQUIRE( phase.name() == "AqueousPhase" );
    REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
    REQUIRE( phase.species().size() == 3 );
    REQUIRE( phase.species(0).name() == "H2O(aq)" );
    REQUIRE( phase.species(1).name() == "H+"      );
    REQUIRE( phase.species(2).name() == "OH-"     );

    auto thermoprops = phase.props(300.00, 10.0e5);

    REQUIRE( thermoprops.temperature() == 300.00 );
    REQUIRE( thermoprops.pressure()    == 10.0e5 );

    REQUIRE( (thermoprops.standardGibbsEnergies()        == 1.0).all() );
    REQUIRE( (thermoprops.standardGibbsEnergies()        == 1.0).all() );
    REQUIRE( (thermoprops.standardEnthalpies()           == 2.0).all() );
    REQUIRE( (thermoprops.standardVolumes()              == 3.0).all() );
    REQUIRE( (thermoprops.standardHeatCapacitiesConstP() == 4.0).all() );
    REQUIRE( (thermoprops.standardHeatCapacitiesConstV() == 5.0).all() );

    ArrayXr n(3); n << 55.508, 1.0e-7, 1.0e-7;

    auto chemicalprops = phase.props(310.00, 11.0e5, n);

    REQUIRE( chemicalprops.temperature() == 310.00 );
    REQUIRE( chemicalprops.pressure()    == 11.0e5 );

    REQUIRE( (chemicalprops.lnActivities()           == 123.0).all() );
    REQUIRE( (chemicalprops.lnActivityCoefficients() == 123.0).all() );
}
