import numpy as np
import os
import pytest
import pandas as pd
import sys

from collections import namedtuple

#pytest
from pytest_regressions.plugin import num_regression

#PyReaktoro
from reaktoro import * 

#some tools
from kineticProblemsSetup import *
sys.path.append(os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir)))
from pythonTools import *


timePropetie = namedtuple('timePropetie', ['ti', 'tf', 'unit'])

@pytest.mark.parametrize('setup, timePropetie, checkedVariables',
    [
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite'),
            timePropetie(0, 24, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite'),
            timePropetie(0, 48, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite'),
            timePropetie(0, 72, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_HCl_CaCO3_Calcite'),
            timePropetie(0, 5, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_HCl_CaCO3_Calcite'),
            timePropetie(0, 10, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
        (
            pytest.lazy_fixture('equilibriumKinectSetUp_H2O_HCl_CaCO3_Calcite'),
            timePropetie(0, 20, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
    ],
    ids = [
        'kinetic problem with H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite 0 to 24 hours',
        'kinetic problem with H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite 0 to 48 hours',
        'kinetic problem with H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite 0 to 72 hours',
        'kinetic problem with H2O_HCl_CaCO3_Calcite 0 to 5 minutes',
        'kinetic problem with H2O_HCl_CaCO3_Calcite 0 to 10 minutes',
        'kinetic problem with H2O_HCl_CaCO3_Calcite 0 to 20 minutes',
        
    ]                                         
)
def test_kinetic_path_solver(num_regression,
                             tmpdir,
                             setup,
                             timePropetie,
                             checkedVariables
                             ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a kinetic problem 
    @param setup
        a tuple that has some objects from kineticProblemSetup.py
        (state, reactions, partition)
    @param timePropetie
        time information about the kinetic problem.
        timePropetie.ti = initial time
        timePropetie.tf = final time
        timePropetie.unit = ti and tf units
    @param checkedVariables 
        a list that has all the variables that will be tested
    '''
    (state, reactions, partition) = setup 
    
    path = KineticPath(reactions)
    
    path.setPartition(partition)
     
    output = path.output()
    output.filename(tmpdir.dirname+"/kinetictPathResult.txt")
    for checkedVariable in checkedVariables:
        output.add(checkedVariable)     
         
    path.solve(state, timePropetie.ti, timePropetie.tf, timePropetie.unit)
     
    pathKineticTable = pd.read_csv(tmpdir.dirname+"/kinetictPathResult.txt", index_col=None, skiprows=1, delim_whitespace=True)
    pathKineticTable.columns = checkedVariables
     
    pathKinectDic = TableToDictionary(pathKineticTable) 
    
    num_regression.check(pathKinectDic)