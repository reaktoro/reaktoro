import numpy as np
import os
import pytest
import pandas as pd
import sys

#pytest
from pytest_regressions.plugin import num_regression

#PyReaktoro
from equilibriumProblemsSetup import *
from PyReaktoro import *

sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir)))
from pythonTools import *


def test_EquilibriumPathSolve(
    num_regression,
    tmpdir,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium path between two states   
    '''
    
    database = Database(b'supcrt98.xml')
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase(b'H O C Na Cl')
    
    system = ChemicalSystem(editor)
    
    problem1 = EquilibriumProblem(system)
    problem1.add(b'H2O', 1, b'kg');
    problem1.add(b'CO2', 0.5, b'mol');
    problem1.add(b'HCl', 1, b'mol');
    
    problem2 = EquilibriumProblem(system)
    problem2.add(b'H2O', 1, b'kg');
    problem2.add(b'CO2', 0.5, b'mol');
    problem2.add(b'NaOH', 2, b'mol');
        
    state1 = equilibrate(problem1)
    state2 = equilibrate(problem2)
    
    path = EquilibriumPath(system)

    output = path.output()
    output.filename(tmpdir.dirname+'/equilibriumPathResult.txt')

    #Define which outputs will be written and checked
    output.add(b't')
    output.add(b'pH')
    output.add(b'speciesMolality(HCO3-)')
    output.add(b'speciesMolality(CO2(aq))')
    output.add(b'speciesMolality(CO3--)')
    
    path.solve(state1, state2)
        
    pathTable = pd.read_csv(tmpdir.dirname+"/equilibriumPathResult.txt", index_col=None, delim_whitespace=True)
    
    pathDict = TableToDictionary(pathTable) 
    
    num_regression.check(pathDict)