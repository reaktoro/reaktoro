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

from reaktoro import ChemicalEditor
import pytest

def _CheckChemicalSystems(system_1, system_2):
    assert system_1.numElements() == system_2.numElements()
    assert system_1.numSpecies() == system_2.numSpecies()
    assert system_1.numPhases() == system_2.numPhases()

def test_add_phases_right_use():
    """Test the normal use of addAqueousPhase, addGaseousPhase and addMineralPhase."""
    editor_first_style = ChemicalEditor()
    editor_first_style.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--")
    editor_first_style.addGaseousPhase("H2O(g) CO2(g)")
    editor_first_style.addMineralPhase("Graphite")
    
    editor_second_style = ChemicalEditor()
    editor_second_style.addAqueousPhase(["H2O(l)", "H+", "OH-", "HCO3-", "CO2(aq)", "CO3--"])
    editor_second_style.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor_second_style.addMineralPhase(["Graphite"])
 
    _CheckChemicalSystems(editor_first_style.createChemicalSystem(), editor_second_style.createChemicalSystem())
    
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
    editor_first_style = ChemicalEditor()
    editor_first_style.addAqueousPhaseWithElements("H C O Ca")
    editor_first_style.addGaseousPhaseWithElements("H C O")
    editor_first_style.addMineralPhaseWithElements("Ca C O")
    
    editor_second_style = ChemicalEditor()
    editor_second_style.addAqueousPhaseWithElements(["H", "C", "O", "Ca"])
    editor_second_style.addGaseousPhaseWithElements(["H", "C", "O"])
    editor_second_style.addMineralPhaseWithElements(["Ca", "C", "O"])
    
    _CheckChemicalSystems(editor_first_style.createChemicalSystem(), editor_second_style.createChemicalSystem())
    
def test_add_phases_with_elements_wrong_use():
    """Test the wrong usage of addAqueousPhaseWithElements."""
    editor = ChemicalEditor()
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements("H2O(l) C Ca")
        
    with pytest.raises(RuntimeError):
        editor.addAqueousPhaseWithElements(["H2O C Ca"])
    
def test_add_phases_with_elements_of_right_use():
    """Test the normal use of addAqueousPhaseWithElementsOf, addGaseousPhaseWithElementsOf and addMineralPhaseWithElementsOf."""
    editor_first_style = ChemicalEditor()
    editor_first_style.addAqueousPhaseWithElementsOf("H2O Ca")
    editor_first_style.addGaseousPhaseWithElementsOf("CO2 H")
    editor_first_style.addMineralPhaseWithElementsOf("CaCO3")
    
    editor_second_style = ChemicalEditor()
    editor_second_style.addAqueousPhaseWithElementsOf(["H2O", "Ca"])
    editor_second_style.addGaseousPhaseWithElementsOf(["CO2", "H"])
    editor_second_style.addMineralPhaseWithElementsOf(["CaCO3"])
    
    _CheckChemicalSystems(editor_first_style.createChemicalSystem(), editor_second_style.createChemicalSystem())
    
    
    