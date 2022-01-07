# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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

    system = ChemicalSystem(phases)

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

    state.setSpeciesAmount(0, 3.1)
    assert state.speciesAmount(0) == 3.1

    state.setSpeciesAmount(0, 3.1, "mmol")
    assert state.speciesAmount(0)[0] == pytest.approx(0.0031)

    state.setSpeciesAmount(N - 1, 4.3)
    assert state.speciesAmount(N - 1) == 4.3

    state.setSpeciesAmount(N - 1, 4.3, "mmol")
    assert state.speciesAmount(N - 1)[0] == pytest.approx(0.0043)

    state.setSpeciesAmount("H2O(aq)", 5.0)
    assert state.speciesAmount("H2O(aq)") == 5.0

    state.setSpeciesAmount("H2O(aq)", 5.0, "mmol")
    assert state.speciesAmount("H2O(aq)")[0] == pytest.approx(0.005)

    state.setSpeciesAmount("Quartz", 10.0)
    assert state.speciesAmount("Quartz") == 10.0

    state.setSpeciesAmount("Quartz", 10.0, "mmol")
    assert state.speciesAmount("Quartz")[0] == pytest.approx(0.01)

    state.setSpeciesMass(0, 7.1)
    assert state.speciesMass(0) == 7.1

    state.setSpeciesMass(0, 7.1, "mg")
    assert state.speciesMass(0)[0] == pytest.approx(7.1e-6)

    state.setSpeciesMass(N - 1, 9.3)
    assert state.speciesMass(N - 1) == 9.3

    state.setSpeciesMass(N - 1, 9.3, "mg")
    assert state.speciesMass(N - 1)[0] == pytest.approx(9.3e-6)

    state.setSpeciesMass("H2O(aq)", 5.0)
    assert state.speciesMass("H2O(aq)") == 5.0

    state.setSpeciesMass("H2O(aq)", 5.0, "mg")
    assert state.speciesMass("H2O(aq)")[0] == pytest.approx(5e-6)

    state.setSpeciesMass("Quartz", 10.0)
    assert state.speciesMass("Quartz") == 10.0

    state.setSpeciesMass("Quartz", 10.0, "mg")
    assert state.speciesMass("Quartz")[0] == pytest.approx(10e-6)

    state.setSpeciesAmount("Calcite", 0.0)
    state.add("Calcite", 1.0, "mol")
    assert state.speciesAmount("Calcite") == 1.0

    state.add("Calcite", 2.0, "mol")
    assert state.speciesAmount("Calcite") == 3.0

    state.add("Calcite", 5.0, "mmol")
    assert state.speciesAmount("Calcite") == 3.005

    state.setSpeciesAmount("Calcite", 0.0)
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

    props = state.props()
    props.update(state)

    assert props.temperature() == state.temperature()
    assert props.pressure() == state.pressure()
    assert props.speciesAmounts() == state.speciesAmounts()
    assert state.charge() == -3.0
