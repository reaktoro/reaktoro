import numpy as np
import os
import pandas as pd
import pytest
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir)))
from pythonTools import StateToDictionary, TableToDictionary
from reaktoro import ChemicalEditor, ChemicalSystem, Database, equilibrate, EquilibriumPath, EquilibriumProblem 


def test_EquilibriumPathSolve(
    num_regression,
    tmpdir,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium path between two states   
    '''
    
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase("H O C Na Cl")
    
    system = ChemicalSystem(editor)
    
    problem1 = EquilibriumProblem(system)
    problem1.add("H2O", 1, "kg");
    problem1.add("CO2", 0.5, "mol");
    problem1.add("HCl", 1, "mol");
    
    problem2 = EquilibriumProblem(system)
    problem2.add("H2O", 1, "kg");
    problem2.add("CO2", 0.5, "mol");
    problem2.add("NaOH", 2, "mol");
        
    state1 = equilibrate(problem1)
    state2 = equilibrate(problem2)
    
    path = EquilibriumPath(system)

    output = path.output()
    output.filename(tmpdir.dirname+"/equilibriumPathResult.txt")

    #Define which outputs will be written and checked
    output.add("t")
    output.add("pH")
    output.add("speciesMolality(HCO3-)")
    output.add("speciesMolality(CO2(aq))")
    output.add("speciesMolality(CO3--)")
    
    path.solve(state1, state2)
        
    pathTable = pd.read_csv(tmpdir.dirname+"/equilibriumPathResult.txt", index_col=None, delim_whitespace=True)
    
    pathDict = TableToDictionary(pathTable) 
    
    num_regression.check(pathDict)