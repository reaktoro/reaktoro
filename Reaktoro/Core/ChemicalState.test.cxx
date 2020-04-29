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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ChemicalState class", "[ChemicalState]")
{
    ActivityPropsFn activity_props_fn = [](ActivityPropsRef props, ActivityArgs args) {};

    // Create the ChemicalSystem object
    ChemicalSystem system({
        Phase()
            .withName("AqueousSolution")
            .withSpecies(SpeciesList("H2O(aq) H+(aq) OH-(aq) H2(aq) O2(aq) Na+(aq) Cl-(aq) NaCl(aq) HCO3-(aq) CO2(aq) CO3--(aq)"))
            .withStateOfMatter(StateOfMatter::Liquid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("GaseousSolution")
            .withSpecies(SpeciesList("H2O(g) CO2(g) H2(g) O2(g)"))
            .withStateOfMatter(StateOfMatter::Gas)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("Halite")
            .withSpecies(SpeciesList("NaCl(s)"))
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("CaCO3(s)")
            .withSpecies(SpeciesList("CaCO3(s)"))
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("Quartz")
            .withSpecies(SpeciesList("SiO2(s)"))
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn)
    });

    ChemicalState state(system);

    auto idx = [&](String name) { return system.species().index(name); };

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setTemperature
    //-------------------------------------------------------------------------
    state.setTemperature(300.0);
    CHECK( state.temperature() == 300.0 );

    state.setTemperature(30.0, "celsius");
    CHECK( state.temperature() == 273.15 + 30.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setPressure
    //-------------------------------------------------------------------------
    state.setPressure(123.0e5);
    CHECK( state.pressure() == 123.0e5 );

    state.setPressure(30.0, "bar");
    CHECK( state.pressure() == 30.0 * 1e5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmounts(value)
    //-------------------------------------------------------------------------
    state.setSpeciesAmounts(1.0);
    CHECK( (state.speciesAmounts() == 1.0).all() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmounts(values)
    //-------------------------------------------------------------------------
    state.setSpeciesAmounts( ArrayXr::Ones(system.species().size()) * 11.0 );
    CHECK( (state.speciesAmounts() == 11.0).all() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(ispecies, amount)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount(3, 4.0);
    CHECK( state.speciesAmount(3) == 4.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(ispecies, amount, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount(3, 4.0, "mmol");
    CHECK( state.speciesAmount(3) == 0.004 );
    CHECK( state.speciesAmount(3, "mmol") == 4.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(name, amount)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount("CO2(g)", 7.0);
    CHECK( state.speciesAmount(idx("CO2(g)")) == 7.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(name, amount, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount("CaCO3(s)", 9.0, "kmol");
    CHECK( state.speciesAmount(idx("CaCO3(s)")) == 9000.0 );
    CHECK( state.speciesAmount(idx("CaCO3(s)"), "kmol") == 9.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(ispecies, mass)
    //-------------------------------------------------------------------------
    state.setSpeciesMass(5, 10.0);
    CHECK( state.speciesMass(5) == 10.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(ispecies, mass, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesMass(5, 15.0, "g");
    CHECK( state.speciesMass(5) == 0.015 );
    CHECK( state.speciesMass(5, "g") == 15.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(name, mass)
    //-------------------------------------------------------------------------
    state.setSpeciesMass("CO2(g)", 7.0);
    CHECK( state.speciesMass(idx("CO2(g)")) == 7.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(name, mass, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesMass("CaCO3(s)", 7.0, "g");
    CHECK( state.speciesMass(idx("CaCO3(s)")) == 0.007 );
    CHECK( state.speciesMass(idx("CaCO3(s)"), "g") == 7.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::phaseProps(iphase)
    //-------------------------------------------------------------------------
    state.setTemperature(400.0);
    state.setPressure(123.0e5);
    state.setSpeciesAmounts(1.0);

    for(auto i = 0; i < system.phases().size(); ++i)
    {
        auto phaseprops = state.phaseProps(i);
        REQUIRE(  phaseprops.temperature()    == 400.0      );
        REQUIRE(  phaseprops.pressure()       == 123.0e5    );
        REQUIRE( (phaseprops.speciesAmounts() == 1.0).all() );
    }

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::props()
    //-------------------------------------------------------------------------
    state.setTemperature(500.0);
    state.setPressure(453.0e5);
    state.setSpeciesAmounts(3.0);

    auto props = state.props();
    REQUIRE(  props.temperature()    == 500.0      );
    REQUIRE(  props.pressure()       == 453.0e5    );
    REQUIRE( (props.speciesAmounts() == 3.0).all() );
}
