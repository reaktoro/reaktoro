import pytest
import numpy as np
import sys
import os

#pytest
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

#PyReaktoro
from equilibriumProblemsSetup import *
from reaktoro import *

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from pythonTools import *

#DO NOT TRY TO PUT ALL EQUILIBRIUM TEST IN A SINGUE TEST
#if one one of then fail the test will stop and you won't test all
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
def test_Equilibrate_problem(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using 
    equilibrate(const EquilibriumProblem& problem) 
    and
    equilibrate(const EquilibriumInverseProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    
    (system, problem) = setup
    
    equilibriumState = ChemicalState()
    
    equilibriumState = equilibrate(problem)
    
    stateDict = StateToDictionary(equilibriumState)
    
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )
def test_Equilibritrate_problem_options(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using 
    equilibrate(const EquilibriumProblem& problem, const EquilibriumOptions& options) 
    and
    equilibrate(const EquilibriumInverseProblem& problem, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    options = EquilibriumOptions()
    
    equilibriumState = equilibrate(problem, options)
     
    stateDict = StateToDictionary(equilibriumState)
    
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3 and Magnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )
def test_Equilibrite_state(
        setup,
        num_regression,
        ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using 
    equilibrate(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    
    (system, problem) = setup
    
    #find a state based on the equilibrium setup
    equilibriumState = equilibrate(problem)
    
    #change pressure and temperature
    equilibriumState.setTemperature(problem.temperature()+20)
    equilibriumState.setPressure(problem.pressure()+20)
    
    #compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(equilibriumState)
     
    stateDict = StateToDictionary(equilibriumState)
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )    
def test_Equilibrate_state_partition(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using 
    equilibrate(ChemicalState& state, const Partition& partition)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem) 
    '''
    (system, problem) = setup
    
    #find a state based on the equilibrium setup  
    state = equilibrate(problem)
    
    #change pressure and temperature
    state.setTemperature(problem.temperature()+10)
    state.setPressure(problem.pressure()+20)
     
    partition = Partition(system)
    
    #compute equilibrium for state with new temperature and pressure 
    equilibriumResult = equilibrate(state, partition)
       
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2, NaCl and Halite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )    
def test_Equilibrate_state_option(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using 
    equilibrate(ChemicalState& state, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    #Find a state based on the equilibrium setup  
    state = equilibrate(problem)
    
    #Change pressure and temperature
    state.setTemperature(problem.temperature()+10)
    state.setPressure(problem.pressure()+20)
     
    options = EquilibriumOptions()

    #Compute equilibrium for state with new temperature and pressure 
    equilibriumResult = equilibrate(state, options)
     
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )
def test_Equilibrate_state_partition_option(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using     
    equilibrate(ChemicalState& state, const Partition& partition, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    #find a state based on the equilibrium setup  
    state = equilibrate(problem)
    
    #change pressure and temperature
    state.setTemperature(problem.temperature()+10)
    state.setPressure(problem.pressure()+20)
     
    options = EquilibriumOptions()
    partition = Partition(system)
     
    #compute equilibrium for state with new temperature and pressure 
    equilibriumResult = equilibrate(state, partition, options)
     
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )
def test_Equilibrate_state_problem(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using       
    equilibrate(ChemicalState& state, const EquilibriumProblem& problem)
    and
    equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = ChemicalState(system)
    
    equilibriumResult = equilibrate(state, problem)    
    
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
    ids=['Eq Problem-H2O,CO2,NaClandHalite at 60C and 300 bar',
        'Eq Problem-H2O,CO2,NaClandHalite already dissolved at 60C and 300 bar',
        'Eq Problem-H2O,Fe(OH)2,Fe(OH)3,NH3andMagnetite',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with fixed amount, activity and pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,Mg,C with defined pH',
        'Eq Inv Problem-H,O,Na,Cl,Ca,C,Calcite with defined pH and fixed amount',
        'Eq Inv Problem-H2O,CO2,NaCl,CaCO3,Calcite with fixed species mass and amounts',
        'Eq Inv Problem - fixed mass, amount and alkalinity',
        'Eq Inv Problem - fixed phase volume'
         ]
    )
def test_Equilibrate_state_problem_option(
    setup,
    num_regression,
    ):
    '''
    An integration test that checks result's reproducibility of 
    the calculation of a problem using      
    equilibrate(ChemicalState& state, const EquilibriumProblem& problem, const EquilibriumOptions& options)
    and
    equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    '''
    (system, problem) = setup
    
    state = ChemicalState(system)
     
    options = EquilibriumOptions()
 
    equilibriumResult = equilibrate(state, problem, options)
     
    stateDict = StateToDictionary(state)
      
    num_regression.check(stateDict,
                         default_tolerance=dict(atol=1e-5, rtol=1e-16))