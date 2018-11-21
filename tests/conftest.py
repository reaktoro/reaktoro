import pytest

from python_tools import convert_table_to_dict, convert_reaktoro_state_to_dict
from reaktoro import Database, ChemicalEditor, ChemicalSystem, EquilibriumProblem, EquilibriumInverseProblem 

@pytest.fixture(scope='function')
def equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar():
    '''
    Build a problem with 1 kg of H2O, 100 g of CO2 and 0.1 mol of NaCl 
    at 60 °C and 300 bar 
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


@pytest.fixture(scope='function')
def equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar():
    '''
    Build a problem with H2O, H+, Na+, Cl-, HCO3-, CO2(aq), CO3-- and 
    Halite at 60 °C and 300 bar 
    '''
    database = Database("supcrt98.xml")
    
    editor = ChemicalEditor(database)
    aqueous = editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--", "CO(aq)"])
    aqueous.setActivityModelDrummondCO2()
    gaseous = editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    gaseous.setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase("Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 100, "g")
    problem.add("NaCl", 1, "mol")
    problem.setTemperature(60, "celsius")
    problem.setPressure(300, "bar")
    
    return (system, problem)

@pytest.fixture(scope='function')
def equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite():
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
    
@pytest.fixture(scope='function')
def equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity():
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

@pytest.fixture(scope='function')
def equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph():
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

@pytest.fixture(scope='function')
def equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts():
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
    
@pytest.fixture(scope='function')
def equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass():
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

@pytest.fixture(scope='function')
def equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity():
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

@pytest.fixture(scope='function')
def equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume():
    '''
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite 
    with fixed values of phase volume 
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

@pytest.fixture(scope='function')
def state_regression(num_regression):
    class StateRegression:
        def check(self, state, tol=None, default_tol=None):
            num_regression.check(
                                 convert_reaktoro_state_to_dict(state), 
                                 basename=None,
                                 fullpath=None, 
                                 tolerances=tol, 
                                 default_tolerance=default_tol, 
                                 data_index=None, 
                                 fill_different_shape_with_nan=True
                                 )
            
    return StateRegression()

@pytest.fixture(scope='function')
def table_regression(num_regression):
    class TableRegression:
        def check(self, table, tol=None, default_tol=None):
            num_regression.check(
                                 convert_table_to_dict(table), 
                                 basename=None,
                                 fullpath=None, 
                                 tolerances=tol, 
                                 default_tolerance=default_tol, 
                                 data_index=None, 
                                 fill_different_shape_with_nan=True
                                 )
                        
    return TableRegression()
