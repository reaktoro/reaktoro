import pytest
from reaktoro import *

"""
All the tests tries to ensure the capability of Reaktoro on solving system that
has liquidlike phases and single mixture
"""
from _pytest.python_api import approx

@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (153.15, 10132.5   , [1.0,0.0]),
        (153.15, 101325.0  , [1.0,0.0]),
        (153.15, 2026500.0 , [0.0,1.0]),
        (153.15, 3039750.0 , [0.0,1.0]),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all CH4 should be gas",
        "temperature equal 153.15 K and 101325 Pa - all CH4 should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all CH4 should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all CH4 should be liquid",
    ],
)
def test_liquid_like_with_CH4(temperature, pressure, result):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was CH4 which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(g)
        result[1] = amount of CH4(oil)
    """
    
    db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    """
    Current version of Reaktoro can only work with system that has H2O, it was
    added with amount equal 0
    """ 
    editor.addAqueousPhase(["H2O(l)"])
    
    phases = []
    
    phases.append(editor.convertAqueousPhase(editor.aqueousPhase()))
    
    #Add gas phase
    gas_species = []
    gas_species.append(db.gaseousSpecies("CH4(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = HydrocarbonPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertHydrocarbonPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("CH4(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, 
                             problem, 
                             phaseIdentificationMethod.Gibbs_residual_based, 
                             0)
    
    assert (state.speciesAmount("CH4(g)"), 
            state.speciesAmount("CH4(oil)")) == approx(result)



@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (323.15, 1013250.0 , [1.0,0.0]),
        (323.15, 2026500.0 , [1.0,0.0]),
        (323.15, 5066250.0 , [0.0,1.0]),
        (323.15, 6079500.0 , [0.0,1.0]),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 101325 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all H2S should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all H2S should be liquid",
    ],
)
def test_liquid_like_with_H2S(temperature, pressure, result):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of H2S(g)
        result[1] = amount of H2S(oil)
    """
    
    db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    """
    Current version of Reaktoro can only work with system that has H2O, it was
    added with amount equal 0
    """ 
    editor.addAqueousPhase(["H2O(l)"])
    
    phases = []
    
    phases.append(editor.convertAqueousPhase(editor.aqueousPhase()))
    
    #Add gas phase
    gas_species = []
    gas_species.append(db.gaseousSpecies("H2S(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = HydrocarbonPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertHydrocarbonPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("H2S(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, 
                             problem, 
                             phaseIdentificationMethod.Gibbs_residual_based, 
                             0)
    
    assert (state.speciesAmount("H2S(g)"), 
            state.speciesAmount("H2S(oil)")) == approx(result) 
    
    

@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (273.15, 1000000.0  , [1.0,0.0]),
        (273.15, 2000000.0  , [1.0,0.0]),
        (273.15, 4000000.0 , [0.0,1.0]),
        (273.15, 5000000.0 , [0.0,1.0]),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 101325 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all H2S should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all H2S should be liquid",
    ],
)
def test_liquid_like_with_CO2(temperature, pressure, result):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CO2(g)
        result[1] = amount of CO2(oil)
    """
    
    db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    """
    Current version of Reaktoro can only work with system that has H2O, it was
    added with amount equal 0
    """ 
    editor.addAqueousPhase(["H2O(l)"])
    
    phases = []
    
    phases.append(editor.convertAqueousPhase(editor.aqueousPhase()))
    
    #Add gas phase
    gas_species = []
    gas_species.append(db.gaseousSpecies("CO2(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = HydrocarbonPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertHydrocarbonPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("CO2(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, 
                             problem, 
                             phaseIdentificationMethod.Gibbs_residual_based, 
                             0)
    
    assert (state.speciesAmount("CO2(g)"), 
            state.speciesAmount("CO2(oil)")) == approx(result)