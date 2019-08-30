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

import pytest
import re
from reaktoro import (
    ChemicalEditor, 
    ChemicalSystem,
    Database,
    AqueousPhase,
    GaseousPhase,
    LiquidPhase,
    MineralPhase,
    PhaseType
)


def _check_equivalent_chemical_systems(system_1, system_2):
    assert system_1.numElements() == system_2.numElements()
    assert system_1.numSpecies() == system_2.numSpecies()
    assert system_1.numPhases() == system_2.numPhases()

def _getting_species_names(phase):
    species_names = ""
    for specie in phase.species():
        species_names += specie.name() + " "
    
    return species_names

def test_add_phases_right_use():
    """Test the normal use of addAqueousPhase, addGaseousPhase and addMineralPhase."""
    database = Database("supcrt98.xml")
    
    # Check adding phase by giving an string with all species
    list_of_aqueous_species_expected = r"H2O\(l\)\s*H\+\s*OH-\s*HCO3-\s*CO2\(aq\)\s*CO3--\s*"
    list_of_gaseous_species_expected = r"H2O\(g\)\s*CO2\(g\)\s*"
    list_of_liquid_species_expected = r"H2O\(liq\)\s*CO2\(liq\)\s*"
    list_of_mineral_species_expected = r"Graphite\s*"
    
    editor1 = ChemicalEditor(database)
    editor1.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--")
    editor1.addGaseousPhase("H2O(g) CO2(g)")
    editor1.addLiquidPhase("H2O(liq) CO2(liq)")
    editor1.addMineralPhase("Graphite")
    
    aqueous_phase_1 = editor1.aqueousPhase()
    gaseous_phase_1 = editor1.gaseousPhase()
    liquid_phase_1 = editor1.liquidPhase()
    mineral_phases_1 = editor1.mineralPhases()    
    
    aqueous_species_added_1 = _getting_species_names(aqueous_phase_1)
    gaseous_species_added_1 = _getting_species_names(gaseous_phase_1)    
    liquid_species_added_1 = _getting_species_names(liquid_phase_1)
    mineral_species_added_1 = _getting_species_names(mineral_phases_1[0])
        
    assert re.match(list_of_aqueous_species_expected, aqueous_species_added_1)
    assert re.match(list_of_gaseous_species_expected, gaseous_species_added_1)
    assert re.match(list_of_liquid_species_expected, liquid_species_added_1)
    assert re.match(list_of_mineral_species_expected, mineral_species_added_1)
    
    # Check adding phase by giving a list of string with all species
    editor2 = ChemicalEditor(database)
    editor2.addAqueousPhase(["H2O(l)", "H+", "OH-", "HCO3-", "CO2(aq)", "CO3--"])
    editor2.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor2.addLiquidPhase(["H2O(liq)", "CO2(liq)"])
    editor2.addMineralPhase(["Graphite"]) 

    aqueous_phase_2 = editor2.aqueousPhase()
    gaseous_phase_2 = editor2.gaseousPhase()
    liquid_phase_2 = editor2.liquidPhase()
    mineral_phases_2 = editor2.mineralPhases()
    
    aqueous_species_added_2 = _getting_species_names(aqueous_phase_2)
    gaseous_species_added_2 = _getting_species_names(gaseous_phase_2)    
    liquid_species_added_2 = _getting_species_names(liquid_phase_2)
    mineral_species_added_2 = _getting_species_names(mineral_phases_2[0])
    
    assert re.match(list_of_aqueous_species_expected, aqueous_species_added_2)
    assert re.match(list_of_gaseous_species_expected, gaseous_species_added_2)
    assert re.match(list_of_liquid_species_expected, liquid_species_added_2)
    assert re.match(list_of_mineral_species_expected, mineral_species_added_2)
    
def test_add_phases_wrong_use():
    """Test the wrong usage of addAqueousPhase, addGaseousPhase and addMineralPhase."""
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    with pytest.raises(RuntimeError):
        editor.addAqueousPhase("H2O(l) C Ca")
        
    with pytest.raises(RuntimeError):
        editor.addAqueousPhase(["H2O C Ca"])
        
    with pytest.raises(RuntimeError):
        editor.addGaseousPhase("CO2(g) H")
        
    with pytest.raises(RuntimeError):
        editor.addMineralPhase("Siderite C")
        
    with pytest.raises(RuntimeError):
        editor.addMineralPhase(["CaCO3"])
        
def test_add_phases_with_elements_right_use():
    """Test the normal use of addAqueousPhaseWithElements, addGaseousPhaseWithElements and addMineralPhaseWithElements."""
    database = Database("supcrt98.xml")
    
    editor1 = ChemicalEditor(database)
    editor1.addAqueousPhaseWithElements("H C O Ca")
    editor1.addGaseousPhaseWithElements("H C O")
    
    editor1.addMineralPhaseWithElements("Ca C O")
    
    editor2 = ChemicalEditor(database)
    editor2.addAqueousPhaseWithElements(["H", "C", "O", "Ca"])
    editor2.addGaseousPhaseWithElements(["H", "C", "O"])
    editor2.addMineralPhaseWithElements(["Ca", "C", "O"])
    
    system1 = ChemicalSystem(editor1)
    system2 = ChemicalSystem(editor2)
    
    _check_equivalent_chemical_systems(system1, system2)
    
def test_add_phases_with_elements_wrong_use():
    """Test the wrong usage of addAqueousPhaseWithElements."""
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements("H2O(l) C Ca")
        
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements(["H2O C Ca"])
    
def test_add_phases_with_elements_of_right_use():
    """Test the normal use of addAqueousPhaseWithElementsOf, addGaseousPhaseWithElementsOf and addMineralPhaseWithElementsOf."""
    database = Database("supcrt98.xml")
    
    editor1 = ChemicalEditor(database)
    editor1.addAqueousPhaseWithElementsOf("H2O Ca")
    editor1.addGaseousPhaseWithElementsOf("CO2 H")
    editor1.addMineralPhaseWithElementsOf("CaCO3")
    
    editor2 = ChemicalEditor(database)
    editor2.addAqueousPhaseWithElementsOf(["H2O", "Ca"])
    editor2.addGaseousPhaseWithElementsOf(["CO2", "H"])
    editor2.addMineralPhaseWithElementsOf(["CaCO3"])
    
    system1 = ChemicalSystem(editor1)
    system2 = ChemicalSystem(editor2)
    
    _check_equivalent_chemical_systems(system1, system2)

def test_chemical_editor_adding_and_getting_phases():
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    gaseous_phase = GaseousPhase()
    liquid_phase = LiquidPhase()
    mineral_phase = MineralPhase()
    
    editor.addPhase(gaseous_phase)
    editor.addPhase(liquid_phase)
    editor.addPhase(mineral_phase)
    
    assert editor.gaseousPhase().name() == gaseous_phase.name()
    assert editor.liquidPhase().name() == liquid_phase.name()
    assert editor.mineralPhases()[0].name() == mineral_phase.name()

def test_chemical_editor_create_system():
    expected = ["H2O(l)", "H+", "OH-", "H2O(g)",
                 "CO2(g)", "H2O(liq)", "CO2(liq)",
                 "Graphite" ]
        
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase("H2O(l) H+ OH-")
    editor.addGaseousPhase("H2O(g) CO2(g)")
    editor.addLiquidPhase("H2O(liq) CO2(liq)")
    editor.addMineralPhase("Graphite")    
    
    system = ChemicalSystem(editor)
    
    species_name = []
    for specie in system.species():
        species_name.append(specie.name())
    
    assert species_name == expected
    
    
