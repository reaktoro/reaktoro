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
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O NaCl CO2")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Halite")
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    
    return (system, problem)


@pytest.fixture(scope='session')
def equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P():
    '''
    Build a problem with H2O, H+, Na+, Cl-, HCO3-, CO2(aq), CO3-- and 
    Halite at 60ºC and 300 bar 
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase([b"H2O(l)", b"H+", b"OH-", b"Na+", b"Cl-", b"HCO3-", b"CO2(aq)", b"CO3--", b"CO(aq)"]) \
        .setActivityModelDrummondCO2()
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"]). \
        setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase(b"Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 1, b"mol")
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite():
    '''
    Build a problem with H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O Fe(OH)2 Fe(OH)3 NH3")
    editor.addGaseousPhase(b"NH3(g)")
    editor.addMineralPhase(b"Magnetite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"Fe(OH)2", 1, b"mol")
    problem.add(b"Fe(OH)3", 2, b"mol")
    problem.add(b"NH3", 1, b"mmol")
    
    return (system, problem)
    
@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with fixed
    species amount, activity and defined pH  
    '''
    
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H O Na Cl Ca Mg C")
    editor.addGaseousPhase(b"H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.add(b"CaCl2", 2, b"mmol")
    problem.add(b"MgCl2", 4, b"mmol")
    problem.pH(3.0, b"HCl")
    problem.fixSpeciesAmount(b"CO2(g)", 1.0, b"mol")
    problem.fixSpeciesActivity(b"O2(g)", 0.20)
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CpH():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with defined pH  
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H O Na Cl Ca Mg C")
    editor.addGaseousPhase(b"H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.add(b"CaCl2", 2, b"mmol")
    problem.add(b"MgCl2", 4, b"mmol")
    problem.pH(4.0, b"CO2")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount():
    '''
    Build a problem with H, O, Na, Cl, Ca, C and Calcite with defined pH 
    and fixed species amount  
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H O Na Cl Ca C")
    editor.addMineralPhase(b"Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.pH(8.0, b"HCl", b"NaOH")
    problem.fixSpeciesAmount(b"Calcite", 1, b"mol")

    return (system, problem)
    
@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupH2O_NaCl_CaCO3_CalcilteFixedMass():
    '''
    Build a problem with H2O, NaCL, CaCO3, CO2, Calcite with fixed
    species mass and amount  
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O NaCl CaCO3")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.fixSpeciesMass(b"Calcite", 100, b"g")
    problem.fixSpeciesAmount(b"CO2(g)", 1.0, b"mol")
    
    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumInverseProblemSetupFixedMAssAmountAndAlkalinity():
    '''
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite 
    with fixed values of Species Mass, Amount and alkalinity 
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H2O NaCl CaCO3")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Calcite")
  
    system = ChemicalSystem(editor)
  
    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.fixSpeciesMass(b"Calcite", 100, b"g")
    problem.fixSpeciesAmount(b"CO2(g)", 1.0, b"mol")
    problem.alkalinity(25.0, b"meq/L", b"Cl")

    return (system, problem)

@pytest.fixture(scope='session')
def equilibriumIverseProblemSetupFixedPhaseVolume():
    '''
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite 
    with fixed values of Phase volume 
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H2O NaCl CaCO3")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Calcite")
  
    system = ChemicalSystem(editor)
  
    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.fixPhaseVolume(b"Gaseous", 0.2, b"m3", b"CO2")
    problem.fixPhaseVolume(b"Aqueous", 0.3, b"m3", b"1 kg H2O; 0.1 mol NaCl")
    problem.fixPhaseVolume(b"Calcite", 0.5, b"m3", b"CaCO3")
    
    return (system, problem)