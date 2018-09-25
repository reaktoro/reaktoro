import pytest
import numpy as np
import sys
import os

# pytest
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

# PyReaktoro
from equilibriumProblemsSetup import *
from PyReaktoro import *

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from pythonTools import *


@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CpH')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH2O_NaCl_CaCO3_CalcilteFixedMass')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupFixedMAssAmountAndAlkalinity')),
        (pytest.lazy_fixture('equilibriumIverseProblemSetupFixedPhaseVolume'))
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite',
        'Eq Inv Problem - H, O, Na, Cl, Ca, Mg, C with fixed amount, activity and pH',
        'Eq Inv Problem - H, O, Na, Cl, Ca, Mg, C with defined pH',
        'Eq Inv Problem - H, O, Na, Cl, Ca, C, Calcite with defined pH and fixed amount',
        'Eq Inv Problem - H2O, CO2, NaCl, CaCO3, Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
    ]
    )
def test_Solve_state(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::solve(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = equilibrate(problem)
    state.setTemperature(problem.temperature() + 30)
    state.setPressure(problem.pressure() + 40)
     
    solver = EquilibriumSolver(system)

    solver.solve(state)

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite'
        ]
    )
def test_Solve_state_problem(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state, problem)

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite'
        ]
    )
def test_Solve_state_T_P_Be(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state,
                problem.temperature(),
                problem.pressure(),
                problem.elementAmounts())

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))
    

@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_Mg_CpH')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupH2O_NaCl_CaCO3_CalcilteFixedMass')),
        (pytest.lazy_fixture('equilibriumInverseProblemSetupFixedMAssAmountAndAlkalinity')),
        (pytest.lazy_fixture('equilibriumIverseProblemSetupFixedPhaseVolume'))
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite',
        'Eq Inv Problem - H, O, Na, Cl, Ca, Mg, C with fixed amount, activity and pH',
        'Eq Inv Problem - H, O, Na, Cl, Ca, Mg, C with defined pH',
        'Eq Inv Problem - H, O, Na, Cl, Ca, C, Calcite with defined pH and fixed amount',
        'Eq Inv Problem - H2O, CO2, NaCl, CaCO3, Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
    ]
    )
def test_Approximate_state(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::approximate(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = equilibrate(problem)
    state.setTemperature(problem.temperature() + 30)
    state.setPressure(problem.pressure() + 40)

    solver = EquilibriumSolver(system)

    solver.approximate(state)

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite'
        ]
    )
def test_Approximate_state_problem(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::approximate(ChemicalState& state, const EquilibriumProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state, problem)

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('setup',
    [
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')),
        (pytest.lazy_fixture('equilibriumProblemSetupH2O_FeOH2_FeOH3_NH3_Magnetite')),
    ],
    ids=[
        'Eq Problem - H2O, CO2, NaCl and Halite at 60C and 300 bar',
        'Eq Problem - H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem - H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite'
        ]
    )
def test_Approximate_state_T_P_be(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of an equilibrium of a state using 
    EquilibriumSolver::approximate(ChemicalState& state, double T, double P, VectorConstRef be)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state,
                        problem.temperature(),
                        problem.pressure(),
                        problem.elementAmounts())

    stateDict = StateToDictionary(state)
    
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))
