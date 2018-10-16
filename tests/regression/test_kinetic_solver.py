import numpy as np
import pandas as pd
import pytest

from collections import namedtuple
from python_tools import convert_dataframe_to_dict
from reaktoro import ChemicalEditor, ChemicalSystem, Database, equilibrate, EquilibriumProblem, KineticPath, Partition, ReactionSystem 


@pytest.fixture(scope='session')
def kinect_problem_with_h2o_hcl_caco3_calcite():
    
    database = Database("supcrt98.xml")
        
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O HCl CaCO3")
    editor.addMineralPhase("Calcite")
    
    calciteReaction = editor.addMineralReaction("Calcite")
    calciteReaction.setEquation("Calcite = Ca++ + CO3--")
    calciteReaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
    calciteReaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
    calciteReaction.setSpecificSurfaceArea(10, "cm2/g")
    
    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)
    
    partition = Partition(system)
    partition.setKineticPhases(["Calcite"])
    
    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add("H2O", 1, "kg")
    problem.add("HCl", 1, "mmol")
    
    state = equilibrate(problem)
    
    state.setSpeciesMass("Calcite", 100, "g")

    return (state, reactions, partition)

@pytest.fixture(scope='session')
def kinect_problem_with_h2o_nacl_caco3_mgco3_calcite_magnesite_dolomite_halite():
        
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase("H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Magnesite")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Halite")
    
    calciteReaction = editor.addMineralReaction("Calcite")
    calciteReaction.setEquation("Calcite = Ca++ + CO3--")
    calciteReaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
    calciteReaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
    calciteReaction.setSpecificSurfaceArea(10, "cm2/g")
        
    magnesiteReaction.addMineralReaction("Magnesite")
    magnesiteReaction.setEquation("Magnesite = Mg++ + CO3--")
    magnesiteReaction.addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol")
    magnesiteReaction.addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
    magnesiteReaction.setSpecificSurfaceArea(10, "cm2/g")
        
    dolomiteReaction.addMineralReaction("Dolomite")
    dolomiteReaction.setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--")
    dolomiteReaction.addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
    dolomiteReaction.addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
    dolomiteReaction.setSpecificSurfaceArea(10, "cm2/g")
        
    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)
    
    partition = Partition(system)
    partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])
    
    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 1, "mol")
    problem.add("CO2", 1, "mol")
    
    
    state = equilibrate(problem)
    
    state.setSpeciesMass("Calcite", 100, "g")
    state.setSpeciesMass("Dolomite", 50, "g")
    
    return (state, reactions, partition)

timePropetie = namedtuple('timePropetie', ['ti', 'tf', 'unit'])

@pytest.mark.xfail
@pytest.mark.parametrize('setup, timePropetie, checkedVariables',
    [
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_nacl_caco3_mgco3_calcite_magnesite_dolomite_halite'),
            timePropetie(0, 24, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_nacl_caco3_mgco3_calcite_magnesite_dolomite_halite'),
            timePropetie(0, 48, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_nacl_caco3_mgco3_calcite_magnesite_dolomite_halite'),
            timePropetie(0, 72, 'hours'),
            ["time(units=hour)","pH","elementMolality(Ca units=molal)", "elementMolality(Mg units=molal)","phaseMass(Calcite units=grams)","phaseMass(Dolomite units=grams)"]
        ),
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_hcl_caco3_calcite'),
            timePropetie(0, 5, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_hcl_caco3_calcite'),
            timePropetie(0, 10, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
        (
            pytest.lazy_fixture('kinect_problem_with_h2o_hcl_caco3_calcite'),
            timePropetie(0, 20, 'minute'),
            ["time(units=minute)", "elementMolality(Ca units=mmolal)", "phaseMass(Calcite units=g)"]
        ),
    ],
    ids=[
        'kinetic prob-h2o nacl caco3 mgco3 co2 calcite magnesite dolomite halite 0 to 24 hours',
        'kinetic prob-h2o nacl caco3 mgco3 co2 calcite magnesite dolomite halite 0 to 48 hours',
        'kinetic prob-h2o nacl caco3 mgco3 co2 calcite magnesite dolomite halite 0 to 72 hours',
        'kinetic prob-h2o hcl caco3 calcite 0 to 5 minutes',
        'kinetic prob-h2o hcl caco3 calcite 0 to 10 minutes',
        'kinetic prob-h2o hcl caco3 calcite 0 to 20 minutes'
        ]                                         
    )
def test_kinetic_path_solve(num_regression,
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
     
    pathKinectDic = convert_dataframe_to_dict(pathKineticTable) 
    
    num_regression.check(pathKinectDic)