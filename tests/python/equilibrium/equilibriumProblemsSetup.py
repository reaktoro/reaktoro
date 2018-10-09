import pytest

from reaktoro import *

#This file has some problems setup that can be used inside pytest
#all fixtures return a tuple with the following structure: (system, problem)

@pytest.fixture(scope='session')
def equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P():
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O NaCl CO2")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Halite")
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 100, "g")
    problem.add("NaCl", 0.1, "mol")
    problem.setTemperature(60, "celsius")
    problem.setPressure(300, "bar")
    
    return (system, problem)


@pytest.fixture(scope='session')
def equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P():
    '''
    Build a problem with H2O, H+, Na+, Cl-, HCO3-, CO2(aq), CO3-- and 
    Halite at 60ºC and 300 bar 
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--", "CO(aq)"]) \
        .setActivityModelDrummondCO2()
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"]). \
        setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase("Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 100, "g")
    problem.add("NaCl", 1, "mol")
    problem.setTemperature(60, "celsius")
    problem.setPressure(300, "bar")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite():
    '''
    Build a problem with H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O Fe(OH)2 Fe(OH)3 NH3")
    editor.addGaseousPhase("NH3(g)")
    editor.addMineralPhase("Magnetite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("Fe(OH)2", 1, "mol")
    problem.add("Fe(OH)3", 2, "mol")
    problem.add("NH3", 1, "mol")
    
    return (system, problem)
    
@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with fixed
    species amount, activity and defined pH  
    '''
    
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H O Na Cl Ca Mg C")
    editor.addGaseousPhase("H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.add("CaCl2", 2, "mmol")
    problem.add("MgCl2", 4, "mmol")
    problem.pH(3.0, "HCl")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
    problem.fixSpeciesActivity("O2(g)", 0.20)
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CpH():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with defined pH  
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H O Na Cl Ca Mg C")
    editor.addGaseousPhase("H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.add("CaCl2", 2, "mmol")
    problem.add("MgCl2", 4, "mmol")
    problem.pH(4.0, "CO2")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount():
    '''
    Build a problem with H, O, Na, Cl, Ca, C and Calcite with defined pH 
    and fixed species amount  
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H O Na Cl Ca C")
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.pH(8.0, "HCl", "NaOH")
    problem.fixSpeciesAmount("Calcite", 1, "mol")

    return (system, problem)
    
@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH2O_NaCl_CaCO3_CalcilteFixedMass():
    '''
    Build a problem with H2O, NaCL, CaCO3, CO2, Calcite with fixed
    species mass and amount  
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixSpeciesMass("Calcite", 100, "g")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupFixedMAssAmountAndAlkalinity():
    '''
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite 
    with fixed values of Species Mass, Amount and alkalinity 
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
  
    system = ChemicalSystem(editor)
  
    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixSpeciesMass("Calcite", 100, "g")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
    problem.alkalinity(25.0, "meq/L", "Cl")

    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumIverseProblemSetupFixedPhaseVolume():
    '''
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite 
    with fixed values of Phase volume 
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
  
    system = ChemicalSystem(editor)
  
    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixPhaseVolume("Gaseous", 0.2, "m3", "CO2")
    problem.fixPhaseVolume("Aqueous", 0.3, "m3", "1 kg H2O; 0.1 mol NaCl")
    problem.fixPhaseVolume("Calcite", 0.5, "m3", "CaCO3")
    
    return (system, problem)