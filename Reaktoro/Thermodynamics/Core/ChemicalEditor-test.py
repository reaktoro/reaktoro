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

from reaktoro import ChemicalEditor, ChemicalSystem
import pytest

def _check_equivalent_chemical_systems(system_1, system_2):
    assert system_1.numElements() == system_2.numElements()
    assert system_1.numSpecies() == system_2.numSpecies()
    assert system_1.numPhases() == system_2.numPhases()

def test_add_phases_right_use():
    """Test the normal use of addAqueousPhase, addGaseousPhase and addMineralPhase."""
    editor1 = ChemicalEditor()
    editor1.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--")
    editor1.addGaseousPhase("H2O(g) CO2(g)")
    editor1.addMineralPhase("Graphite")
    
    editor2 = ChemicalEditor()
    editor2.addAqueousPhase(["H2O(l)", "H+", "OH-", "HCO3-", "CO2(aq)", "CO3--"])
    editor2.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor2.addMineralPhase(["Graphite"])
    
    system1 = ChemicalSystem(editor1)
    system2 = ChemicalSystem(editor2)
 
    _check_equivalent_chemical_systems(system1, system2)
    
def test_add_phases_wrong_use():
    """Test the wrong usage of addAqueousPhase, addGaseousPhase and addMineralPhase."""
    editor = ChemicalEditor()
    with pytest.raises(RuntimeError):
        editor.addAqueousPhase("H2O(l) C Ca")
        
    with pytest.raises(RuntimeError):        
        editor.addAqueousPhase(["H2O C Ca"])
        
    with pytest.raises(RuntimeError):
        editor.addGaseousPhase("CO2(g) H")
        
    with pytest.raises(RuntimeError):
        editor.addGaseousPhase(["CO2"])
        
    with pytest.raises(RuntimeError):
        editor.addMineralPhase("Siderita C")
        
    with pytest.raises(RuntimeError):
        editor.addMineralPhase(["CaCO3"])
        
def test_add_phases_with_elements_right_use():
    """Test the normal use of addAqueousPhaseWithElements, addGaseousPhaseWithElements and addMineralPhaseWithElements."""
    editor1 = ChemicalEditor()
    editor1.addAqueousPhaseWithElements("H C O Ca")
    editor1.addGaseousPhaseWithElements("H C O")
    editor1.addMineralPhaseWithElements("Ca C O")
    
    editor2 = ChemicalEditor()
    editor2.addAqueousPhaseWithElements(["H", "C", "O", "Ca"])
    editor2.addGaseousPhaseWithElements(["H", "C", "O"])
    editor2.addMineralPhaseWithElements(["Ca", "C", "O"])
    
    system1 = ChemicalSystem(editor1)
    system2 = ChemicalSystem(editor2)
    
    _check_equivalent_chemical_systems(system1, system2)
    
def test_add_phases_with_elements_wrong_use():
    """Test the wrong usage of addAqueousPhaseWithElements."""
    editor = ChemicalEditor()
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements("H2O(l) C Ca")
        
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements(["H2O C Ca"])
    
def test_add_phases_with_elements_of_right_use():
    """Test the normal use of addAqueousPhaseWithElementsOf, addGaseousPhaseWithElementsOf and addMineralPhaseWithElementsOf."""
    editor1 = ChemicalEditor()
    editor1.addAqueousPhaseWithElementsOf("H2O Ca")
    editor1.addGaseousPhaseWithElementsOf("CO2 H")
    editor1.addMineralPhaseWithElementsOf("CaCO3")
    
    editor2 = ChemicalEditor()
    editor2.addAqueousPhaseWithElementsOf(["H2O", "Ca"])
    editor2.addGaseousPhaseWithElementsOf(["CO2", "H"])
    editor2.addMineralPhaseWithElementsOf(["CaCO3"])
    
    system1 = ChemicalSystem(editor1)
    system2 = ChemicalSystem(editor2)
    
    _check_equivalent_chemical_systems(system1, system2)
    
    
    