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
    // Create the Database object for the ChemicalSystem
    Database db;
    db.addSpecies( Species("H2O(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("H+(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("OH-(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("H2(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("O2(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("Na+(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("Cl-(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("NaCl(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("HCO3-(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("CO2(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("CO3--(aq)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("H2O(g)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("CO2(g)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("H2(g)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("O2(g)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("NaCl(s)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("CaCO3(s)").withStandardGibbsEnergy(0.0) );
    db.addSpecies( Species("SiO2(s)").withStandardGibbsEnergy(0.0) );

    // Create the ActivityPropsFn of the Phase objects for the ChemicalSystem
    ActivityPropsFn activity_props_fn = [](ActivityPropsRef props, ActivityArgs args) {};

    // Create the Phase objects for the ChemicalSystem
    const Vec<Phase> phases =
    {
        Phase()
            .withName("AqueousPhase")
            .withSpecies(db.speciesWithAggregateState(AggregateState::Aqueous))
            .withStateOfMatter(StateOfMatter::Liquid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("GaseousPhase")
            .withSpecies(db.speciesWithAggregateState(AggregateState::Gas))
            .withStateOfMatter(StateOfMatter::Gas)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("Halite")
            .withSpecies({ db.species().get("NaCl(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("CaCO3(s)")
            .withSpecies({ db.species().get("CaCO3(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn),
        Phase()
            .withName("Quartz")
            .withSpecies({ db.species().get("SiO2(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activity_props_fn)
    };


    // Create the ChemicalSystem object
    ChemicalSystem system(db, phases);

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
}
