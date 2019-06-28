import pytest
from reaktoro import *
import numpy as np

"""
All the tests tries to ensure the capability of Reaktoro on solving system that
has liquidlike phases and single mixture
"""

from _pytest.python_api import approx

@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (190.0, 8.0   , [0.461176, 0.0999923, 0.0388239, 0.400008]),
        (196.0, 15.0  , [0.409528, 0.0661928, 0.0904722, 0.433807]),                      
    ],
    ids=[
        "temperature equal 190.0 K and 8.0 bar",
        "temperature equal 196.0 K and 15.0 bar",
    ],
)
def test_equilibrium_with_CH4_CO2(temperature, pressure, result):
    """
    This test checks the capability of solving binary mixture  
    The selected species were CH4 and CO2 which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(g)
        result[1] = amount of CO2(g)
        result[2] = amount of CH4(oil)
        result[3] = amount of CO2(oil)
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
    gas_species.append(db.gaseousSpecies("CO2(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add hydrocarbon phase
    oil_species = [] #HydrocarbonSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = LiquidPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertLiquidPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.5, "mol")
    problem.add("CO2(g)", 0.5, "mol")
    
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
    
    assert  (state.speciesAmount("CH4(g)"), 
             state.speciesAmount("CO2(g)"),
             state.speciesAmount("CH4(oil)"),
             state.speciesAmount("CO2(oil)")) == approx(result, abs=1e-6)

@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (273.15, 30.0  , [0.489136, 0.334536, 0.0108636, 0.165463]),
        (293.15, 70.0  , [0.455735, 0.280682, 0.0442652, 0.219318]),                      
    ],
    ids=[
        "temperature equal 273.15 K and 30.0 bar",
        "temperature equal 293.15 K and 70.0 bar",
    ],
)
def test_equilibrium_with_CH4_H2S(temperature, pressure, result):
    """
    This test checks the capability of solving binary mixture  
    The selected species were CH4 and H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(g)
        result[1] = amount of H2S(g)
        result[2] = amount of CH4(oil)
        result[3] = amount of H2S(oil)
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
    gas_species.append(db.gaseousSpecies("H2S(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = LiquidPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertLiquidPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.5, "mol")
    problem.add("H2S(g)", 0.5, "mol")
    
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
            state.speciesAmount("H2S(g)"),
            state.speciesAmount("CH4(oil)"),
            state.speciesAmount("H2S(oil)")) == approx(result, abs=1e-6)


@pytest.mark.parametrize(
    "temperature, pressure, result",
    [
        (233.15, 20.0  , [0.579744, 0.103813, 0.0300209, 0.0202560, 0.246187, 0.0199791]),
        (273.15, 40.0  , [0.593993, 0.295721, 0.0467587, 0.0060072, 0.054279, 0.0032413]),
        (293.15, 70.0  , [0.600000, 0.350000, 0.0500000, 0.0000000, 0.000000, 0.0000000]),
    ],
    ids=[
        "temperature equal 233.15 K and 20.0 bar",
        "temperature equal 273.15 K and 40.0 bar",
        "temperature equal 293.15 K and 70.0 bar",
    ],
)
def test_equilibrium_with_CH4_CO2_H2S(temperature, pressure, result):
    """
    This test checks the capability of solving a ternary mixture with
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    @param result
        a list that has the amounts of all species as following:
        result[0] = amount of CH4(g)
        result[1] = amount of H2S(g)
        result[2] = amount of CO2(g)
        result[3] = amount of CH4(oil)
        result[4] = amount of H2S(oil)
        result[5] = amount of CO2(oil)
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
    gas_species.append(db.gaseousSpecies("H2S(g)"))
    gas_species.append(db.gaseousSpecies("CO2(g)"))
    
    gas_mixture = GaseousMixture(gas_species)
    gas_phase = GaseousPhase(gas_mixture)
    gas_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertGaseousPhase(gas_phase))
    
    #Add oil phase
    oil_species = [] #OilSpecies_vector()
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CH4(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("H2S(oil)")))
    oil_species.append(HydrocarbonSpecies(db.gaseousSpecies("CO2(oil)")))
    
    oil_mixture = OilMixture(oil_species)
    oil_phase = LiquidPhase(oil_mixture)
    oil_phase.setChemicalModelSoaveRedlichKwong()
    
    phases.append(editor.convertLiquidPhase(oil_phase))
    
    system = ChemicalSystem(phases)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.60, "mol")
    problem.add("H2S(g)", 0.35, "mol")
    problem.add("CO2(g)", 0.05, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, 
                             problem, 
                             phaseIdentificationMethod.Gibbs_residual_based, 
                             0)
    
    
    assert (state.speciesAmount("CH4(g)"),
            state.speciesAmount("H2S(g)"),
            state.speciesAmount("CO2(g)"),
            state.speciesAmount("CH4(oil)"),
            state.speciesAmount("H2S(oil)"),
            state.speciesAmount("CO2(oil)")) == approx(result, abs=1e-6)
            
