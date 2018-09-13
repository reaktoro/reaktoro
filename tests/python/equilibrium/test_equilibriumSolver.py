import pytest
import numpy as np

#pytest
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

#PyReaktoro
from problemsSetup import *
from PyReaktoro import *


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

    num_regression.check(stateDict(state))


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

    num_regression.check(stateDict(state))


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

    num_regression.check(stateDict(state))


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

    num_regression.check(stateDict(state))


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

    num_regression.check(stateDict(state))


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

    num_regression.check(stateDict(state))

def stateDict(state):
    #TODO - add pH, pE, Eh, alk -- no pybind
    properties = state.properties() 
    system = state.system()
    
    phase_masses = state.properties().phaseMasses().val
    phase_volumes = state.properties().phaseVolumes().val
    
    checkOutPut = {
        'number of phases': np.asarray([system.numPhases()]),
        'temperature': np.asarray([state.temperature()]),
        'pressure': np.asarray([state.pressure()]),
        'molar fraction': np.asarray(state.properties().moleFractions().val),
        #'lnActivity Coefficient' : np.asarray(state.properties().lnActivityCoefficients().val),
        'Activity' : np.asarray(state.properties().chemicalPotentials().val),
        'phase moles' : np.asarray(state.properties().phaseAmounts().val),
        'phases masses' : np.asarray(state.properties().phaseMasses().val),
        'phase molar volumes' : np.asarray(state.properties().phaseMolarVolumes().val),
        'phase volumes' : np.asarray(state.properties().phaseVolumes().val),
        'phase volume fraction' : np.asarray(phase_volumes/sum(phase_volumes)),
        'phase densities' : np.asarray(phase_masses/phase_volumes),
        'phase stability indices' : np.asarray(state.phaseStabilityIndices()),
        'species amounts' : np.asarray(state.speciesAmounts())
        }
    
    for i in range(0, system.numElements()):
        checkOutPut[system.element(i).name()+' amount'] = np.asarray([state.elementAmount(i)]) 
    
    return checkOutPut
