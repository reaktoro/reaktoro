# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *
import numpy as npy
import pytest


def testChemicalState():
    db = SupcrtDatabase("supcrtbl")

    phases = Phases(db)
    phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- Ca+2 Mg+2 HCO3- CO3-2 CO2(aq) SiO2(aq)") )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhase("Halite") )
    phases.add( MineralPhase("Calcite") )
    phases.add( MineralPhase("Magnesite") )
    phases.add( MineralPhase("Dolomite") )
    phases.add( MineralPhase("Quartz") )

    reactions = Reactions()
    reactions.add( db.reaction("H2O(aq) = H+ + OH-") )
    reactions.add( db.reaction("CO2(aq) + H2O(aq) = H+ + HCO3-") )
    reactions.add( db.reaction("CO2(g) = CO2(aq)") )
    reactions.add( db.reaction("Halite = Na+ + Cl-") )
    reactions.add( db.reaction("Calcite") )

    system = ChemicalSystem(phases, reactions)

    state = ChemicalState(system)

    N = system.species().size()
    n = npy.arange(N)

    assert state.temperature() == 298.15
    assert state.pressure() == 1e5
    assert npy.max(state.speciesAmounts()) == 1e-16  # default initial amount of all species

    state.temperature(30.0, "celsius")
    assert state.temperature() == pytest.approx(303.15)

    state.pressure(1.0, "MPa")
    assert state.pressure() == pytest.approx(1e6)

    state.temperature(40, "celsius")  # testing with int argument
    assert state.temperature() == pytest.approx(313.15)

    state.pressure(2, "MPa")  # testing with int argument
    assert state.pressure() == pytest.approx(2e6)

    state.setTemperature(30.0, "celsius")
    assert state.temperature() == pytest.approx(303.15)

    state.setPressure(1.0, "MPa")
    assert state.pressure() == pytest.approx(1e6)

    state.setTemperature(40, "celsius")  # testing with int argument
    assert state.temperature() == pytest.approx(313.15)

    state.setPressure(2, "MPa")  # testing with int argument
    assert state.pressure() == pytest.approx(2e6)

    state.setSpeciesAmounts(1.0)
    assert npy.min(state.speciesAmounts()) == 1.0
    assert npy.max(state.speciesAmounts()) == 1.0
    assert state.charge() == 1.0 # +1 -1 +1 -1 +2 +2 -1 -2 = 1

    state.setSpeciesAmounts(n)
    assert npy.all(state.speciesAmounts() == n)
    assert state.charge() == -3.0 # +1 -1*2 +1*3 -1*4 +2*5 +2*6 -1*7 -2*8 = -3

    state.setSpeciesAmount(0, 3.1, "mol")
    assert state.speciesAmount(0) == 3.1

    state.setSpeciesAmount(0, 3.1, "mmol")
    assert state.speciesAmount(0).val() == pytest.approx(0.0031)

    state.setSpeciesAmount(N - 1, 4.3, "mol")
    assert state.speciesAmount(N - 1) == 4.3

    state.setSpeciesAmount(N - 1, 4.3, "mmol")
    assert state.speciesAmount(N - 1).val() == pytest.approx(0.0043)

    state.setSpeciesAmount("H2O(aq)", 5.0, "mol")
    assert state.speciesAmount("H2O(aq)") == 5.0

    state.setSpeciesAmount("H2O(aq)", 5.0, "mmol")
    assert state.speciesAmount("H2O(aq)").val() == pytest.approx(0.005)

    state.setSpeciesAmount("Quartz", 10.0, "mol")
    assert state.speciesAmount("Quartz") == 10.0

    state.setSpeciesAmount("Quartz", 10.0, "mmol")
    assert state.speciesAmount("Quartz").val() == pytest.approx(0.01)

    state.setSpeciesMass(0, 7.1, "kg")
    assert state.speciesMass(0) == 7.1

    state.setSpeciesMass(0, 7.1, "mg")
    assert state.speciesMass(0).val() == pytest.approx(7.1e-6)

    state.setSpeciesMass(N - 1, 9.3, "kg")
    assert state.speciesMass(N - 1) == 9.3

    state.setSpeciesMass(N - 1, 9.3, "mg")
    assert state.speciesMass(N - 1).val() == pytest.approx(9.3e-6)

    state.setSpeciesMass("H2O(aq)", 5.0, "kg")
    assert state.speciesMass("H2O(aq)") == 5.0

    state.setSpeciesMass("H2O(aq)", 5.0, "mg")
    assert state.speciesMass("H2O(aq)").val() == pytest.approx(5e-6)

    state.setSpeciesMass("Quartz", 10.0, "kg")
    assert state.speciesMass("Quartz") == 10.0

    state.setSpeciesMass("Quartz", 10.0, "mg")
    assert state.speciesMass("Quartz").val() == pytest.approx(10e-6)

    state.setSpeciesAmount("Calcite", 0.0, "mol")
    state.add("Calcite", 1.0, "mol")
    assert state.speciesAmount("Calcite") == 1.0

    state.add("Calcite", 2.0, "mol")
    assert state.speciesAmount("Calcite") == 3.0

    state.add("Calcite", 5.0, "mmol")
    assert state.speciesAmount("Calcite") == 3.005

    state.setSpeciesAmount("Calcite", 0.0, "mol")
    state.add("Calcite", 2000.0, "g")
    assert state.speciesMass("Calcite") == 2.0  # in kg

    state.add("Calcite", 3000.0, "g")
    assert state.speciesMass("Calcite") == 5.0  # in kg

    state.add("Calcite", 5000.0, "mg")
    assert state.speciesMass("Calcite") == 5.005  # in kg

    state.set("Calcite", 3.0, "mol")
    assert state.speciesAmount("Calcite") == 3.0  # in kg

    state.set("Calcite", 5000.0, "mg")
    assert state.speciesMass("Calcite") == 0.005  # in kg

    scaled = ChemicalState(state)
    scaled.scaleSpeciesAmounts(2.0)
    assert (scaled.speciesAmounts().asarray()
        == state.speciesAmounts().asarray() * 2.0).all()

    scaled = ChemicalState(state)
    scaled.scaleSpeciesAmountsInPhase("AqueousPhase", 3.0)
    assert (scaled.speciesAmountsInPhase("AqueousPhase").asarray()
        == state.speciesAmountsInPhase("AqueousPhase").asarray() * 3.0).all()

    scaled = ChemicalState(state)
    scaled.scaleAmount(1.0, "mol")
    scaled.props().update(scaled)
    assert scaled.props().amount().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scalePhaseAmount("AqueousPhase", 1.0, "mol")
    scaled.props().update(scaled)
    assert scaled.props().phaseProps("AqueousPhase").amount().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scalePhaseAmount("GaseousPhase", 1.0, "mol")
    scaled.props().update(scaled)
    assert scaled.props().phaseProps("GaseousPhase").amount().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scaleFluidAmount(2.0, "mol")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithFluidState()
    assert sum([scaled.props().phaseProps(i).amount().val() for i in iphases]) == pytest.approx(2.0)

    scaled = ChemicalState(state)
    scaled.scaleSolidAmount(3.0, "mol")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithSolidState()
    assert sum([scaled.props().phaseProps(i).amount().val() for i in iphases]) == pytest.approx(3.0)

    scaled = ChemicalState(state)
    scaled.scaleMass(1.0, "kg")
    scaled.props().update(scaled)
    assert scaled.props().mass().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scalePhaseMass("AqueousPhase", 1.0, "kg")
    scaled.props().update(scaled)
    assert scaled.props().phaseProps("AqueousPhase").mass().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scaleFluidMass(2.0, "kg")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithFluidState()
    assert sum([scaled.props().phaseProps(i).mass() for i in iphases]) == pytest.approx(2.0)

    scaled = ChemicalState(state)
    scaled.scaleSolidMass(3.0, "kg")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithSolidState()
    assert sum([scaled.props().phaseProps(i).mass() for i in iphases]) == pytest.approx(3.0)

    scaled = ChemicalState(state)
    scaled.scaleVolume(1.0, "m3")
    scaled.props().update(scaled)
    assert scaled.props().volume().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scalePhaseVolume("AqueousPhase", 1.0, "m3")
    scaled.props().update(scaled)
    assert scaled.props().phaseProps("AqueousPhase").volume().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scalePhaseVolume("GaseousPhase", 1.0, "m3")
    scaled.props().update(scaled)
    assert scaled.props().phaseProps("GaseousPhase").volume().val() == pytest.approx(1.0)

    scaled = ChemicalState(state)
    scaled.scaleFluidVolume(2.0, "m3")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithFluidState()
    assert sum([scaled.props().phaseProps(i).volume().val() for i in iphases]) == pytest.approx(2.0)

    scaled = ChemicalState(state)
    scaled.scaleSolidVolume(3.0, "m3")
    scaled.props().update(scaled)
    iphases = scaled.props().indicesPhasesWithSolidState()
    assert sum([scaled.props().phaseProps(i).volume().val() for i in iphases]) == pytest.approx(3.0)

    props = state.props()
    props.update(state)

    assert props.temperature() == state.temperature()
    assert props.pressure() == state.pressure()
    assert props.speciesAmounts() == state.speciesAmounts()
    assert state.charge() == -3.0

    #-------------------------------------------------------------------------
    # TESTING METHOD: ChemicalState.update
    #-------------------------------------------------------------------------

    otherstate = ChemicalState(system)

    otherstate.update(state.temperature(), state.pressure(), state.speciesAmounts())

    assert otherstate.temperature() == state.temperature()
    assert otherstate.pressure() == state.pressure()
    assert otherstate.speciesAmounts() == state.speciesAmounts()

    state.props().update(state)

    assert otherstate.props().temperature() == state.temperature()
    assert otherstate.props().pressure() == state.pressure()
    assert otherstate.props().speciesAmounts() == state.speciesAmounts()
    assert otherstate.props().speciesActivitiesLn() == state.props().speciesActivitiesLn()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ChemicalState.update
    #-------------------------------------------------------------------------

    otherstate.updateIdeal(state.temperature(), state.pressure(), state.speciesAmounts())

    assert otherstate.temperature() == state.temperature()
    assert otherstate.pressure() == state.pressure()
    assert otherstate.speciesAmounts() == state.speciesAmounts()

    state.props().updateIdeal(state)

    assert otherstate.props().temperature() == state.temperature()
    assert otherstate.props().pressure() == state.pressure()
    assert otherstate.props().speciesAmounts() == state.speciesAmounts()
    assert otherstate.props().speciesActivitiesLn() == state.props().speciesActivitiesLn()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ChemicalState.assign
    #-------------------------------------------------------------------------
    state = ChemicalState(system)
    otherstate = ChemicalState(system)

    state.setTemperature(321, "K")
    state.setPressure(432, "MPa")
    state.setSpeciesAmounts(1.23)

    otherstate.assign(state)

    assert otherstate.temperature() == state.temperature()
    assert otherstate.pressure() == state.pressure()
    assert otherstate.speciesAmounts() == state.speciesAmounts()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ChemicalState.clone
    #-------------------------------------------------------------------------
    state = ChemicalState(system)
    state.setTemperature(321, "K")

    otherstate = state.clone()

    assert otherstate.temperature() == state.temperature()

    otherstate.setTemperature(1234, "K")

    assert state.temperature() == 321  # ensure state was not changed with changing otherstate above!
