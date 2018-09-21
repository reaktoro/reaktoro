import os
import pytest
import numpy as np
import pandas as pd

#pytest
from pytest_regressions.plugin import num_regression

#PyReaktoro
from problemsSetup import *
from PyReaktoro import *
from pythonTools import *

def test_equilibrium_path_solve(num_regression,tmpdir):
    
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase(b"H O C Na Cl")
    
    system = ChemicalSystem(editor)
    
    problem1 = EquilibriumProblem(system)
    problem1.add(b"H2O", 1, b"kg");
    problem1.add(b"CO2", 0.5, b"mol");
    problem1.add(b"HCl", 1, b"mol");
    
    problem2 = EquilibriumProblem(system)
    problem2.add(b"H2O", 1, b"kg");
    problem2.add(b"CO2", 0.5, b"mol");
    problem2.add(b"NaOH", 2, b"mol");
        
    state1 = equilibrate(problem1)
    state2 = equilibrate(problem2)
    
    path = EquilibriumPath(system)

    output = path.output()
    output.filename(tmpdir.dirname+"/equilibriumPathResult.txt")
    output.add(b"t")
    output.add(b"pH")
    output.add(b"speciesMolality(HCO3-)", b"HCO3- [molal]")
    output.add(b"speciesMolality(CO2(aq))", b"CO2(aq) [molal]")
    output.add(b"speciesMolality(CO3--)", b"CO3-- [molal]")

    path.solve(state1, state2)
    
    pathTable = pd.read_csv(tmpdir.dirname+"/equilibriumPathResult.txt", index_col=None, skiprows=1, delim_whitespace=True)
    pathTable.columns = ["t", "pH", "HCO3- [molal]", "CO2(aq) [molal]", "CO3-- [molal]"]
    
    pathDic = pathTableDictionary(pathTable) 
    
    
    num_regression.check(pathDic)    
