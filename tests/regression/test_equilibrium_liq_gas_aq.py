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
All the tests tries to ensure the capability of Reaktoro on solving system that
has LiquidPhase, GaseousPhase and AqueousPhase with water
"""
from _pytest.python_api import approx

@pytest.mark.parametrize(
    "temperature, pressure",
    [
        (380.35, 75.6  ),
        (380.35, 122.7 ),
        (310.95, 130.0 ),
        (310.95, 62.6  ),
    ],
    ids=[
        "temperature equal 380.35 K and 75.6 bar",
        "temperature equal 380.35 K and 122.7 bar",
        "temperature equal 310.95 K and 130.0 bar",
        "temperature equal 310.95 K and 62.6 bar",
    ],
)
def test_equilibrium_CH4_H2S_CO2_H2O_liq_gas_aq(temperature, pressure, num_regression):
    """
    This test checks the capability of solving a system that has CH4, H2S,
    CO2, H2O with
    @param Temperature
        temperature in Kelvin which will be used to compute equilibrium
    @param Pressure
        pressure in bar which will be used to compute equilibrium
    """
    
    db = Database("supcrt98.xml") 
     
    editor = ChemicalEditor(db)
    
    editor.addAqueousPhase(["CO2(aq)",  "H2S(aq)", "H2O(l)" ])
    editor.addGaseousPhase(["CH4(g)", "CO2(g)", "H2S(g)", "H2O(g)"])
    editor.addLiquidPhase(["CH4(liq)", "CO2(liq)", "H2S(liq)", "H2O(liq)"])
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    
    
    problem.setTemperature(temperature, "K")
    problem.setPressure(pressure, "bar")
    problem.add("H2O(g)", 0.50, "mol")
    problem.add("CO2(g)", 0.05, "mol")
    problem.add("H2S(g)", 0.40, "mol")
    problem.add("CH4(g)", 0.05, "mol")
    
    solver = EquilibriumSolver(problem.system())
    
    options = EquilibriumOptions();
    options.hessian = GibbsHessian.Exact;
    options.nonlinear.max_iterations = 100;
    options.optimum.max_iterations = 200;
    options.optimum.ipnewton.step = StepMode.Conservative;
    options.optimum.tolerance = 1e-17
    solver.setOptions(options)
            
    state = ChemicalState(system)
    
    res_final = solver.solve(state, problem)
    

    species_amount = {
        "CO2(aq)": np.asarray([state.speciesAmount("CO2(g)")]),
        "H2S(aq)": np.asarray([state.speciesAmount("H2S(aq)")]),
        "H2O(l)": np.asarray([state.speciesAmount("H2O(l)")]),
        "CH4(g)": np.asarray([state.speciesAmount("CH4(g)")]),
        "CO2(g)": np.asarray([state.speciesAmount("CO2(g)")]),
        "H2S(g)": np.asarray([state.speciesAmount("H2S(g)")]),
        "H2O(g)": np.asarray([state.speciesAmount("H2O(g)")]),
        "CH4(liq)": np.asarray([state.speciesAmount("CH4(liq)")]),
        "CO2(liq)": np.asarray([state.speciesAmount("CO2(liq)")]),
        "H2S(liq)": np.asarray([state.speciesAmount("H2S(liq)")]),
        "H2O(liq)": np.asarray([state.speciesAmount("H2O(liq)")]),
        }

    num_regression.check(species_amount)
