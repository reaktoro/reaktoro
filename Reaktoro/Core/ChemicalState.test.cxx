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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

namespace test {

/// Return a mock ChemicalSystem object for test reasons.
auto createChemicalSystem() -> ChemicalSystem;

} // namespace test

TEST_CASE("Testing ChemicalState class", "[ChemicalState]")
{
    ChemicalSystem system = test::createChemicalSystem();

    ChemicalState state(system);

    auto idx = [&](String name) { return system.species().index(name); };

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::temperature
    //-------------------------------------------------------------------------
    state.temperature(200.0);
    CHECK( state.temperature() == 200.0 );

    state.temperature(33.0, "celsius");
    CHECK( state.temperature() == 273.15 + 33.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::pressure
    //-------------------------------------------------------------------------
    state.pressure(234.0e5);
    CHECK( state.pressure() == 234.0e5 );

    state.pressure(40.0, "bar");
    CHECK( state.pressure() == 40.0 * 1e5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::add(name, value, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount("SiO2(s)", 0.0, "mol");
    state.add("SiO2(s)", 1.0, "mol");
    CHECK( state.speciesAmount("SiO2(s)") == Approx(1.0) );
    state.add("SiO2(s)", 2.0, "mol");
    CHECK( state.speciesAmount("SiO2(s)") == Approx(3.0) );
    state.add("SiO2(s)", 5.0, "mmol");
    CHECK( state.speciesAmount("SiO2(s)") == Approx(3.005) );

    state.setSpeciesAmount("SiO2(s)", 0.0, "mol");
    state.add("SiO2(s)", 2000.0, "g");
    CHECK( state.speciesMass("SiO2(s)") == Approx(2.0) ); // in kg
    state.add("SiO2(s)", 3000.0, "g");
    CHECK( state.speciesMass("SiO2(s)") == Approx(5.0) ); // in kg
    state.add("SiO2(s)", 5000.0, "mg");
    CHECK( state.speciesMass("SiO2(s)") == Approx(5.005) ); // in kg

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::set(name, value, unit)
    //-------------------------------------------------------------------------
    state.set("NaCl(s)", 1.0, "mol");
    CHECK( state.speciesAmount("NaCl(s)") == Approx(1.0) );

    state.set("NaCl(s)", 3000.0, "g");
    CHECK( state.speciesMass("NaCl(s)") == Approx(3.0) ); // in kg

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
    state.setSpeciesAmount(3, 4.0, "mol");
    CHECK( state.speciesAmount(3) == 4.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(ispecies, amount, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount(3, 4.0, "mmol");
    CHECK( state.speciesAmount(3) == 0.004 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(name, amount)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount("CO2(g)", 7.0, "mol");
    CHECK( state.speciesAmount(idx("CO2(g)")) == 7.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesAmount(name, amount, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesAmount("CaCO3(s)", 9.0, "kmol");
    CHECK( state.speciesAmount(idx("CaCO3(s)")) == 9000.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(ispecies, mass)
    //-------------------------------------------------------------------------
    state.setSpeciesMass(5, 10.0, "kg");
    CHECK( state.speciesMass(5) == 10.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(ispecies, mass, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesMass(5, 15.0, "g");
    CHECK( state.speciesMass(5) == Approx(0.015) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(name, mass)
    //-------------------------------------------------------------------------
    state.setSpeciesMass("CO2(g)", 7.0, "kg");
    CHECK( state.speciesMass(idx("CO2(g)")) == 7.0 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::setSpeciesMass(name, mass, unit)
    //-------------------------------------------------------------------------
    state.setSpeciesMass("CaCO3(s)", 7.0, "g");
    CHECK( state.speciesMass(idx("CaCO3(s)")) == 0.007 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::charge
    //-------------------------------------------------------------------------
    state.setSpeciesAmounts(1e-16);
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Na+(aq)", 1.0, "mol");
    state.set("Cl-(aq)", 1.0, "mol");
    state.set("Ca++(aq)", 1.0, "mol");

    CHECK( state.charge() == Approx(2.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scaleSpeciesAmounts
    //-------------------------------------------------------------------------
    ChemicalState scaled(system);
    scaled = ChemicalState(state);
    scaled.scaleSpeciesAmounts(2.0);
    CHECK( scaled.speciesAmounts().isApprox(state.speciesAmounts() * 2.0) );
    CHECK_THROWS( scaled.scaleSpeciesAmounts(-1.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scaleSpeciesAmountsInPhase
    //-------------------------------------------------------------------------
    scaled = ChemicalState(state);
    scaled.scaleSpeciesAmountsInPhase(0, 3.0);
    scaled.scaleSpeciesAmountsInPhase("GaseousPhase", 4.0);
    CHECK( scaled.speciesAmountsInPhase(0).isApprox(state.speciesAmountsInPhase(0) * 3.0) );
    CHECK( scaled.speciesAmountsInPhase("GaseousPhase").isApprox(state.speciesAmountsInPhase("GaseousPhase") * 4.0) );
    CHECK_THROWS( scaled.scaleSpeciesAmountsInPhase(0, -1.0) );
    CHECK_THROWS( scaled.scaleSpeciesAmountsInPhase("XYZ", -1.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scaleVolume
    //-------------------------------------------------------------------------
    scaled = ChemicalState(state);
    scaled.scaleVolume(1.2, "m3");
    scaled.props().update(scaled);
    CHECK( scaled.props().volume() == Approx(1.2) );
    CHECK_THROWS( scaled.scaleVolume(-1.2, "m3") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scalePhaseVolume
    //-------------------------------------------------------------------------
    scaled = ChemicalState(state);
    scaled.scalePhaseVolume("AqueousPhase", 1.2, "m3");
    scaled.props().update(scaled);
    CHECK( scaled.props().phaseProps("AqueousPhase").volume() == Approx(1.2) );
    CHECK_THROWS( scaled.scalePhaseVolume("AqueousPhase", -1.2, "m3") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scaleMass
    //-------------------------------------------------------------------------
    scaled = ChemicalState(state);
    scaled.scaleMass(1.2, "kg");
    scaled.props().update(scaled);
    CHECK( scaled.props().mass() == Approx(1.2) );
    CHECK_THROWS( scaled.scaleMass(-1.2, "kg") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::scalePhaseMass
    //-------------------------------------------------------------------------
    scaled = ChemicalState(state);
    scaled.scalePhaseMass("AqueousPhase", 1.2, "kg");
    scaled.props().update(scaled);
    CHECK( scaled.props().phaseProps("AqueousPhase").mass() == Approx(1.2) );
    CHECK_THROWS( scaled.scalePhaseMass("AqueousPhase", -1.2, "kg") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::update
    //-------------------------------------------------------------------------

    ChemicalState otherstate(system);

    otherstate.update(state.temperature(), state.pressure(), state.speciesAmounts());

    CHECK( otherstate.temperature() == state.temperature() );
    CHECK( otherstate.pressure() == state.pressure() );
    CHECK( otherstate.speciesAmounts().isApprox(state.speciesAmounts()) );

    state.props().update(state);

    CHECK( otherstate.props().temperature() == state.temperature() );
    CHECK( otherstate.props().pressure() == state.pressure() );
    CHECK( otherstate.props().speciesAmounts().isApprox(state.speciesAmounts()) );
    CHECK( otherstate.props().speciesActivitiesLn().isApprox(state.props().speciesActivitiesLn()) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::update
    //-------------------------------------------------------------------------

    otherstate.updateIdeal(state.temperature(), state.pressure(), state.speciesAmounts());

    CHECK( otherstate.temperature() == state.temperature() );
    CHECK( otherstate.pressure() == state.pressure() );
    CHECK( otherstate.speciesAmounts().isApprox(state.speciesAmounts()) );

    state.props().updateIdeal(state);

    CHECK( otherstate.props().temperature() == state.temperature() );
    CHECK( otherstate.props().pressure() == state.pressure() );
    CHECK( otherstate.props().speciesAmounts().isApprox(state.speciesAmounts()) );
    CHECK( otherstate.props().speciesActivitiesLn().isApprox(state.props().speciesActivitiesLn()) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalState::props
    //-------------------------------------------------------------------------
    state.props().update(state);
    CHECK(  state.props().temperature() == state.temperature() );
    CHECK(  state.props().pressure() == state.pressure() );
    CHECK( (state.props().speciesAmounts() == state.speciesAmounts()).all() );

    ChemicalProps props = state.props();

    CHECK(  props.temperature() == state.temperature() );
    CHECK(  props.pressure() == state.pressure() );
    CHECK( (props.speciesAmounts() == state.speciesAmounts()).all() );

    auto createChemicalProps = [system]()
    {
        ChemicalState state(system);
        state.setTemperature(288.0);
        state.setPressure(1.3e5);
        state.set("H2O(aq)", 1.0, "kg");
        state.set("Na+(aq)", 1.0, "mol");
        state.set("Cl-(aq)", 1.0, "mol");
        state.set("Ca++(aq)", 1.0, "mol");
        state.set("Mg++(aq)", 1.0, "mol");
        state.set("CO3--(aq)", 2.0, "mol");
        state.set("HCO3-(aq)", 1.234, "mol");
        state.props().update(state);

        CHECK(  state.props().temperature() == 288.0 );
        CHECK(  state.props().pressure() == 1.3e5 );
        CHECK(  state.props().charge() == Approx(-1.234) );

        return state.props();
    };

    props = createChemicalProps();

    CHECK(  props.temperature() == 288.0 );
    CHECK(  props.pressure() == 1.3e5 );
    CHECK(  props.charge() == Approx(-1.234) );
}
