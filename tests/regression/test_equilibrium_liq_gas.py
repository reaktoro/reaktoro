import numpy as np
import pytest

from reaktoro import (
    ChemicalEditor,
    ChemicalState,
    ChemicalSystem,
    Database,
    EquilibriumOptions,
    EquilibriumProblem,
    EquilibriumSolver,
    GibbsHessian,
    StepMode,
)


"""
Those test tries to ensure the capability of Reaktoro on solving system that
have liquid like phases without water
"""

@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (153.15, 10132.5   ),
        (153.15, 101325.0  ),
        (153.15, 2026500.0 ),
        (153.15, 3039750.0 ),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all CH4 should be gas",
        "temperature equal 153.15 K and 101325 Pa - all CH4 should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all CH4 should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all CH4 should be liquid",
    ],
)
def test_equilibrium_CH4_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was CH4 which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["CH4(g)"])
    editor.addLiquidPhase(["CH4(liq)"])
        
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("CH4(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions()
    options.hessian = GibbsHessian.Exact
    options.nonlinear.max_iterations = 100
    options.optimum.max_iterations = 200
    options.optimum.ipnewton.step = StepMode.Conservative
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded
    
    species_amounts = {
        "CH4(g)": np.asarray([state.speciesAmount("CH4(g)")]),
        "CH4(liq)": np.asarray([state.speciesAmount("CH4(liq)")]),
        }
    
    num_regression.check(species_amounts)


@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (323.15, 1013250.0 ),
        (323.15, 2026500.0 ),
        (323.15, 5066250.0 ),
        (323.15, 6079500.0 ),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 101325 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all H2S should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all H2S should be liquid",
    ],
)
def test_equilibrium_H2S_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["H2S(g)"])
    editor.addLiquidPhase(["H2S(liq)"])
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("H2S(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions()
    options.hessian = GibbsHessian.Exact
    options.nonlinear.max_iterations = 100
    options.optimum.max_iterations = 200
    options.optimum.ipnewton.step = StepMode.Conservative
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded
    
    species_amounts = {
        "H2S(g)": np.asarray([state.speciesAmount("H2S(g)")]),
        "H2S(liq)": np.asarray([state.speciesAmount("H2S(liq)")]),
        }
    
    num_regression.check(species_amounts)  
    

@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (273.15, 1000000.0 ),
        (273.15, 2000000.0 ),
        (273.15, 4000000.0 ),
        (273.15, 5000000.0 ),                      
    ],
    ids=[
        "temperature equal 153.15 K and 10132.5 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 101325 Pa - all H2S should be gas",
        "temperature equal 153.15 K and 2026500 Pa - all H2S should be liquid",
        "temperature equal 153.15 K and 3039750 Pa - all H2S should be liquid",
    ],
)
def test_equilibrium_CO2_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving a system that has 1 specie. 
    The selected specie was H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in Pa which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["CO2(g)"])
    editor.addLiquidPhase(["CO2(liq)"])
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "Pa")
    problem.add("CO2(g)", 1.0, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions()
    options.hessian = GibbsHessian.Exact
    options.nonlinear.max_iterations = 100
    options.optimum.max_iterations = 200
    options.optimum.ipnewton.step = StepMode.Conservative
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded

    species_amounts = {
        "CO2(g)": np.asarray([state.speciesAmount("CO2(g)")]),
        "CO2(liq)": np.asarray([state.speciesAmount("CO2(liq)")]),
        }
    
    num_regression.check(species_amounts)


@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (190.0, 8.0),
        (196.0, 15.0),                      
    ],
    ids=[
        "temperature equal 190.0 K and 8.0 bar",
        "temperature equal 196.0 K and 15.0 bar",
    ],
)
def test_equilibrium_CH4_CO2_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving binary mixture  
    The selected species were CH4 and CO2 which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["CH4(g)", "CO2(g)"])
    editor.addLiquidPhase(["CH4(liq)", "CO2(liq)"])
       
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.5, "mol")
    problem.add("CO2(g)", 0.5, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions()
    options.hessian = GibbsHessian.Exact
    options.nonlinear.max_iterations = 100
    options.optimum.max_iterations = 200
    options.optimum.ipnewton.step = StepMode.Conservative
    options.optimum.tolerance = 1e-14
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded
    
    species_amounts = {
        "CH4(g)": np.asarray([state.speciesAmount("CH4(g)")]),
        "CO2(g)": np.asarray([state.speciesAmount("CO2(g)")]),
        "CH4(liq)": np.asarray([state.speciesAmount("CH4(liq)")]),
        "CO2(liq)": np.asarray([state.speciesAmount("CO2(liq)")])
        }
    
    num_regression.check(species_amounts)

@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (273.15, 30.0),
        (293.15, 70.0),                      
    ],
    ids=[
        "temperature equal 273.15 K and 30.0 bar",
        "temperature equal 293.15 K and 70.0 bar",
    ],
)
def test_equilibrium_CH4_H2S_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving binary mixture  
    The selected species were CH4 and H2S which 
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["CH4(g)", "H2S(g)"])
    editor.addLiquidPhase(["CH4(liq)", "H2S(liq)"])
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.5, "mol")
    problem.add("H2S(g)", 0.5, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions()
    options.hessian = GibbsHessian.Exact
    options.nonlinear.max_iterations = 100
    options.optimum.max_iterations = 200
    options.optimum.ipnewton.step = StepMode.Conservative
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded
    
    species_amount = {
        "CH4(g)": np.asarray([state.speciesAmount("CH4(g)")]),
        "H2S(g)": np.asarray([state.speciesAmount("H2S(g)")]),
        "CH4(liq)": np.asarray([state.speciesAmount("CH4(liq)")]),
        "H2S(liq)": np.asarray([state.speciesAmount("H2S(liq)")]),
        }
    
    num_regression.check(species_amount)


@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (233.15, 20.0),
        (273.15, 40.0),
        (293.15, 70.0),
    ],
    ids=[
        "temperature equal 233.15 K and 20.0 bar",
        "temperature equal 273.15 K and 40.0 bar",
        "temperature equal 293.15 K and 70.0 bar",
    ],
)
def test_equilibrium_CH4_CO2_H2S_liq_gas(temperature, pressure, num_regression):
    """
    This test checks the capability of solving a ternary mixture with
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addGaseousPhase(["CH4(g)", "H2S(g)", "CO2(g)"])
    editor.addLiquidPhase(["CH4(liq)", "H2S(liq)", "CO2(liq)"])
        
    system = ChemicalSystem(editor)
    
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
    
    result = solver.solve(state, problem)
    
    assert result.optimum.succeeded
    
    species_amount = {
        "CH4(g)": np.asarray([state.speciesAmount("CH4(g)")]),
        "H2S(g)": np.asarray([state.speciesAmount("H2S(g)")]),
        "CO2(g)": np.asarray([state.speciesAmount("CO2(g)")]),
        "CH4(liq)": np.asarray([state.speciesAmount("CH4(liq)")]),
        "H2S(liq)": np.asarray([state.speciesAmount("H2S(liq)")]),
        "CO2(liq)": np.asarray([state.speciesAmount("CO2(liq)")]),
        }
    
    num_regression.check(species_amount)
