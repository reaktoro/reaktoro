import pytest

from reaktoro import *

@pytest.fixture(scope='session')
def equilibriumKinectSetUp_H2O_HCl_CaCO3_Calcite():
    
    database = Database("supcrt98.xml")
        
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O HCl CaCO3")
    editor.addMineralPhase("Calcite")
    
    editor.addMineralReaction("Calcite") \
        .setEquation("Calcite = Ca++ + CO3--") \
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, "cm2/g")
    
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
def equilibriumKinectSetUp_H2O_NaCl_CaCO3_MgCO3_CO2_Calcite_Magnesite_Dolomite_Halite():
        
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    
    editor.addAqueousPhase("H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Magnesite")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Halite")
    
    editor.addMineralReaction("Calcite") \
        .setEquation("Calcite = Ca++ + CO3--") \
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, "cm2/g")
        
    editor.addMineralReaction("Magnesite") \
        .setEquation("Magnesite = Mg++ + CO3--") \
        .addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, "cm2/g")
        
    editor.addMineralReaction("Dolomite") \
        .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
        .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
        .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
        .setSpecificSurfaceArea(10, "cm2/g")
        
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