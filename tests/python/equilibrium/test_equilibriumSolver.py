import pytest
import numpy as np

#pytest
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

#PyReaktoro
from problemsSetup import *
from PyReaktoro import *
from pythonTools import *


#DO NOT TRY TO PUT ALL EQUILIBRIUM TEST IN A SINGUE TEST
#if one one of then fail the test will stop and you won't test all
@pytest.mark.parametrize('problemSetup',
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
def test_solve_state_T_P_be(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    system = problemSetup.system()

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state,
                problemSetup.temperature(),
                problemSetup.pressure(),
                problemSetup.elementAmounts())

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))
    


@pytest.mark.parametrize('problemSetup',
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
def test_solve_state_problem(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    system = problemSetup.system()

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state, problemSetup)

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('problemSetup',
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
def test_solve_state(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    
    system = problemSetup.system()
    
    state = equilibrate(problemSetup)
    state.setTemperature(problemSetup.temperature()+30)
    state.setPressure(problemSetup.pressure()+40)

    
    solver = EquilibriumSolver(system)

    solver.solve(state)

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('problemSetup',
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
def test_approximate_state_T_P_be(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    system = problemSetup.system()
    
    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state,
                        problemSetup.temperature(),
                        problemSetup.pressure(),
                        problemSetup.elementAmounts())

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('problemSetup',
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
def test_approximate_state_problem(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    system = problemSetup.system()
    
    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state, problemSetup)

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize('problemSetup',
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
def test_approximate_state(
        problemSetup,
        num_regression):
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    system = problemSetup.system()
    
    state = equilibrate(problemSetup)
    state.setTemperature(problemSetup.temperature() + 30)
    state.setPressure(problemSetup.pressure() + 40)

    solver = EquilibriumSolver(system)

    solver.approximate(state)

    stateDic = stateDictionary(state)
    
    num_regression.check(stateDic, 
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))

