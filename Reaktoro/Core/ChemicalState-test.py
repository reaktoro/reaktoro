# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
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
from pytest import approx, raises
from numpy import array

def test_ChemicalState():

    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--".split())
    editor.addGaseousPhase("H2O(g) CO2(g)".split())
    editor.addMineralPhase("Graphite")
    system = ChemicalSystem(editor)
    state = ChemicalState(system)

    # A sensible array of species amounts
    n = array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0])

    # Set temperature to 300 K
    state.setTemperature(300)
    assert state.temperature() == 300

    # Set temperature to 30 degC (= 303.15 K)
    state.setTemperature(30, "celsius")
    assert state.temperature() == 303.15

    # Set temperature to 330 K
    state.setTemperature(330, "kelvin")
    assert state.temperature() == 330


    # Set pressure to 1e6 Pa
    state.setPressure(1e6)
    assert state.pressure() == 1e6

    # Set pressure to 2e6 Pa
    state.setPressure(2e6, "Pa")
    assert state.pressure() == 2e6

    # Set pressure to 1 kPa (= 1e3 Pa)
    state.setPressure(1, "kPa")
    assert state.pressure() == 1e3

    # Set pressure to 1 MPa (= 1e6 Pa)
    state.setPressure(1, "MPa")
    assert state.pressure() == 1e6

    # Set pressure to 100 bar (= 100e5 Pa)
    state.setPressure(100, "bar")
    assert state.pressure() == 100e5


    # Set all species amounts to 1 mol
    state.setSpeciesAmounts(1.0)
    assert all(state.speciesAmounts() == 1.0)

    # Set all species amounts to a given list of values
    state.setSpeciesAmounts(n)
    assert state.speciesAmounts() == approx(n)

    # Set the amounts of species HCO3-, CO2(aq), CO3-- using their indices
    state.setSpeciesAmounts([0.2, 0.4, 0.02], [3, 4, 5])
    assert state.speciesAmount("HCO3-") == approx(0.2)
    assert state.speciesAmount("CO2(aq)") == approx(0.4)
    assert state.speciesAmount("CO3--") == approx(0.02)

    # Set the amounts of species H+ to 1e-8 mol using its index
    state.setSpeciesAmount(1, 1e-8)
    assert state.speciesAmount("H+") == approx(1e-8)

    # Set the amounts of species H+ to 2e-8 mol using its name
    state.setSpeciesAmount("H+", 2e-8)
    assert state.speciesAmount("H+") == approx(2e-8)

    # Set the amounts of species HCO3- to 1 mmol (= 1e-3 mol) using its index
    state.setSpeciesAmount(3, 1.0, "mmol")
    assert state.speciesAmount("HCO3-") == approx(1e-3)

    # Set the amounts of species HCO3- to 100 umol (= 100e-6 mol) using its name
    state.setSpeciesAmount("HCO3-", 100.0, "umol")
    assert state.speciesAmount("HCO3-") == approx(100e-6)

    # Set the amounts of species H+ to 2e-8 mol using its name
    state.setSpeciesAmount("H+", 2e-8)
    assert state.speciesAmount("H+") == approx(2e-8)

    # Assert a runtime error is thrown when list of amounts has mismatched dimensions
    with raises(RuntimeError):
        state.setSpeciesAmounts([])  # wrong dimension of array


    # Set the mass of species H2O(l) to 1 kg using its index
    state.setSpeciesMass(0, 1.0)
    assert state.speciesAmount(0) * system.species(0).molarMass() == approx(1.0)

    # Set the mass of species H2O(l) to 2 kg using its name
    state.setSpeciesMass("H2O(l)", 2.0)
    assert state.speciesAmount(0) * system.species(0).molarMass() == approx(2.0)

    # Set the mass of species H2O(l) to 1000 g (= 1 kg) using its index
    state.setSpeciesMass(0, 1000, "g")
    assert state.speciesAmount(0) * system.species(0).molarMass() == approx(1.0)

    # Set the mass of species H2O(l) to 2000 g (= 2 kg) using its name
    state.setSpeciesMass("H2O(l)", 2000, "g")
    assert state.speciesAmount(0) * system.species(0).molarMass() == approx(2.0)


    # Set the amounts of all species to 1 mol
    state.setSpeciesAmounts(1.0)

    # Check the amounts in mol of each element using their indices
    assert state.elementAmount(0) == approx(5)   # C
    assert state.elementAmount(1) == approx(7)   # H
    assert state.elementAmount(2) == approx(13)  # O
    assert state.elementAmount(3) == approx(-3)  # Z

    # Check the amounts in mol of each element using their symbols
    assert state.elementAmount("C") == approx(5)
    assert state.elementAmount("H") == approx(7)
    assert state.elementAmount("O") == approx(13)
    assert state.elementAmount("Z") == approx(-3)

    # Check the amounts in mmol of each element using their indices
    assert state.elementAmount(0, "mmol") == approx(5e3)   # C
    assert state.elementAmount(1, "mmol") == approx(7e3)   # H
    assert state.elementAmount(2, "mmol") == approx(13e3)  # O
    assert state.elementAmount(3, "mmol") == approx(-3e3)  # Z

    # Check the amounts in umol of each element using their symbols
    assert state.elementAmount("C", "umol") == approx(5e6)
    assert state.elementAmount("H", "umol") == approx(7e6)
    assert state.elementAmount("O", "umol") == approx(13e6)
    assert state.elementAmount("Z", "umol") == approx(-3e6)


    # Check the amounts in mol of each element in the aqueous phase using element indices
    assert state.elementAmountInPhase(0, 0) == approx(3)   # C
    assert state.elementAmountInPhase(1, 0) == approx(5)   # H
    assert state.elementAmountInPhase(2, 0) == approx(10)  # O
    assert state.elementAmountInPhase(3, 0) == approx(-3)  # Z

    # Check the amounts in mol of each element in the aqueous phase using element symbols
    assert state.elementAmountInPhase("C", "Aqueous") == approx(3)
    assert state.elementAmountInPhase("H", "Aqueous") == approx(5)
    assert state.elementAmountInPhase("O", "Aqueous") == approx(10)
    assert state.elementAmountInPhase("Z", "Aqueous") == approx(-3)

    # Check the amounts in mmol of each element in the aqueous phase using element indices
    assert state.elementAmountInPhase(0, 0, "mmol") == approx(3e3)   # C
    assert state.elementAmountInPhase(1, 0, "mmol") == approx(5e3)   # H
    assert state.elementAmountInPhase(2, 0, "mmol") == approx(10e3)  # O
    assert state.elementAmountInPhase(3, 0, "mmol") == approx(-3e3)  # Z

    # Check the amounts in umol of each element in the aqueous phase using element symbols
    assert state.elementAmountInPhase("C", "Aqueous", "umol") == approx(3e6)
    assert state.elementAmountInPhase("H", "Aqueous", "umol") == approx(5e6)
    assert state.elementAmountInPhase("O", "Aqueous", "umol") == approx(10e6)
    assert state.elementAmountInPhase("Z", "Aqueous", "umol") == approx(-3e6)

    # Check the amounts in mol of each element in the gaseous phase using element indices
    assert state.elementAmountInPhase(0, 1) == approx(1)   # C
    assert state.elementAmountInPhase(1, 1) == approx(2)   # H
    assert state.elementAmountInPhase(2, 1) == approx(3)   # O
    assert state.elementAmountInPhase(3, 1) == approx(0)   # Z

    # Check the amounts in mol of each element in the gaseous phase using element symbols
    assert state.elementAmountInPhase("C", "Gaseous") == approx(1)
    assert state.elementAmountInPhase("H", "Gaseous") == approx(2)
    assert state.elementAmountInPhase("O", "Gaseous") == approx(3)
    assert state.elementAmountInPhase("Z", "Gaseous") == approx(0)

    # Check the amounts in mmol of each element in the gaseous phase using element indices
    assert state.elementAmountInPhase(0, 1, "mmol") == approx(1e3)   # C
    assert state.elementAmountInPhase(1, 1, "mmol") == approx(2e3)   # H
    assert state.elementAmountInPhase(2, 1, "mmol") == approx(3e3)   # O
    assert state.elementAmountInPhase(3, 1, "mmol") == approx(0)     # Z

    # Check the amounts in umol of each element in the gaseous phase using element symbols
    assert state.elementAmountInPhase("C", "Gaseous", "umol") == approx(1e6)
    assert state.elementAmountInPhase("H", "Gaseous", "umol") == approx(2e6)
    assert state.elementAmountInPhase("O", "Gaseous", "umol") == approx(3e6)
    assert state.elementAmountInPhase("Z", "Gaseous", "umol") == approx(0)


    # Check the amount of each phase using their indices
    assert state.phaseAmount(0) == approx(6)
    assert state.phaseAmount(1) == approx(2)

    # Check the amount of each phase using their names
    assert state.phaseAmount("Aqueous") == approx(6)
    assert state.phaseAmount("Gaseous") == approx(2)

    # Check the amount of each phase in mmol using their indices
    assert state.phaseAmount(0, "mmol") == approx(6e3)
    assert state.phaseAmount(1, "mmol") == approx(2e3)

    # Check the amount of each phase in umol using their names
    assert state.phaseAmount("Aqueous", "umol") == approx(6e6)
    assert state.phaseAmount("Gaseous", "umol") == approx(2e6)

    # Set species amounts to sensible values
    state.setSpeciesAmounts(n)

    # Calculate the chemical properties of the state
    properties = state.properties()

    # Check the usage state.properties().someProperty() works
    assert state.properties().phaseVolumes().val == approx(properties.phaseVolumes().val)
    assert state.properties().lnActivities().val == approx(properties.lnActivities().val)
    assert state.properties().chemicalPotentials().val == approx(properties.chemicalPotentials().val)

    # Get the volumes of the phases
    volumes = properties.phaseVolumes().val

    # Scale the species amounts by a factor of 2
    state.scaleSpeciesAmounts(2.0)
    assert state.speciesAmounts() == approx(2 * n)

    # Scale the species amounts in each phase by a factor of 0.5
    state.scaleSpeciesAmountsInPhase(0, 0.5)
    state.scaleSpeciesAmountsInPhase(1, 0.5)
    state.scaleSpeciesAmountsInPhase(2, 0.5)
    assert state.speciesAmounts() == approx(n)

    # Scale the species amounts in each phase using their indices such that each phase have 1 m3
    state.scalePhaseVolume(0, 1.0)
    state.scalePhaseVolume(1, 1.0)
    state.scalePhaseVolume(2, 1.0)
    assert state.properties().phaseVolumes().val == approx(1.0)

    # Scale the species amounts in each phase using their indices such that each phase have 1000 ml (= 1e-3 m3)
    state.scalePhaseVolume(0, 1000.0, "ml")
    state.scalePhaseVolume(1, 1000.0, "ml")
    state.scalePhaseVolume(2, 1000.0, "ml")
    assert state.properties().phaseVolumes().val == approx(1.0e-3)

    # Scale the species amounts in each phase using their names such that each phase have 1 m3
    state.scalePhaseVolume("Aqueous", 1.0)
    state.scalePhaseVolume("Gaseous", 1.0)
    state.scalePhaseVolume("Graphite", 1.0)
    assert state.properties().phaseVolumes().val == approx(1.0)

    # Scale the species amounts in each phase using their names such that each phase have 1000 ml (= 1e-3 m3)
    state.scalePhaseVolume("Aqueous", 1000.0, "ml")
    state.scalePhaseVolume("Gaseous", 1000.0, "ml")
    state.scalePhaseVolume("Graphite", 1000.0, "ml")
    assert state.properties().phaseVolumes().val == approx(1.0e-3)

    # Scale the species amounts in each fluid phase such that the total fluid volume is 1 m3
    state.scaleFluidVolume(1.0)
    assert state.properties().fluidVolume().val == approx(1.0)

    # Scale the species amounts in each fluid phase such that the total fluid volume is 1 cm3 (= 1e-6 m3)
    state.scaleFluidVolume(1.0, "cm3")
    assert state.properties().fluidVolume().val == approx(1.0e-6)

    # Scale the species amounts in each solid phase such that the total solid volume is 1 m3
    state.scaleSolidVolume(1.0)
    assert state.properties().solidVolume().val == approx(1.0)

    # Scale the species amounts in each solid phase such that the total solid volume is 1 cm3 (= 1e-6 m3)
    state.scaleSolidVolume(100.0, "cm3")
    assert state.properties().solidVolume().val == approx(100.0e-6)

    # Scale the species amounts such that the total system volume is 1 cm3 (= 1e-6 m3)
    state.scaleVolume(1.0)
    assert state.properties().volume().val == approx(1.0)

    # Scale the species amounts such that the total system volume is 100 ml (= 100e-6 m3)
    state.scaleVolume(100.0, "ml")
    assert state.properties().volume().val == approx(100.0e-6)
