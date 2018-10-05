import pytest

from reaktoro import *

@pytest.fixture(scope='session')
def equilibriumKinectSetUp_H2O_HCl_CaCO3_Calcite():
    
    database = Database(b"supcrt98.xml")
        
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O HCl CaCO3")
    editor.addMineralPhase(b"Calcite")
    
    editor.addMineralReaction(b"Calcite") \
        .setEquation(b"Calcite = Ca++ + CO3--") \
        .addMechanism(b"logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")
    
    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)
    
    partition = Partition(system)
    partition.setKineticPhases([b"Calcite"])
    
    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"HCl", 1, b"mmol")
    
    state = equilibrate(problem)
    
    state.setSpeciesMass(b"Calcite", 100, b"g")

    return (state, reactions, partition)

@pytest.fixture(scope='session')
def equilibriumKinectSetUp_H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite():
        
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase(b"H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Calcite")
    editor.addMineralPhase(b"Magnesite")
    editor.addMineralPhase(b"Dolomite")
    editor.addMineralPhase(b"Halite")
    
    editor.addMineralReaction(b"Calcite") \
        .setEquation(b"Calcite = Ca++ + CO3--") \
        .addMechanism(b"logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")
        
    editor.addMineralReaction(b"Magnesite") \
        .setEquation(b"Magnesite = Mg++ + CO3--") \
        .addMechanism(b"logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")
        
    editor.addMineralReaction(b"Dolomite") \
        .setEquation(b"Dolomite = Ca++ + Mg++ + 2*CO3--") \
        .addMechanism(b"logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
        .addMechanism(b"logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
        .setSpecificSurfaceArea(10, b"cm2/g")
        
    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)
    
    partition = Partition(system)
    partition.setKineticSpecies([b"Calcite", b"Magnesite", b"Dolomite"])
    
    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 1, b"mol")
    problem.add(b"CO2", 1, b"mol")
    
    
    state = equilibrate(problem)
    
    state.setSpeciesMass(b"Calcite", 100, b"g")
    state.setSpeciesMass(b"Dolomite", 50, b"g")
    
    return (state, reactions, partition)