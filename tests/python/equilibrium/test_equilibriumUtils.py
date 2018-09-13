import pytest
import numpy as np

#pytest
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

#PyReaktoro
from problemsSetup import *
from PyReaktoro import *

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


      
# #TODO: move problem setup to equilibriumProblems after fix xfail
#@pytest.mark.skip(reason="RES-9, some times is doesn't converge, throwing an error")
def test_demo_equilibrium_fixed_alkalinity(num_regression):
    '''
    Test a problem with H2O, NaCl, CO2, CaCO3 and Calcite 
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
  
    state = equilibrate(problem)
  
    output = stateDict(state)
      
    num_regression.check(output, 
                         default_tolerance=dict(atol=1e-7, rtol=1e-18))
  
  
#TODO: add documentation
#TODO: move problem setup to equilibriumProblems after fix xfail
#@pytest.mark.skip(reason="RES-10 some times it doesn't converge")
def test_demo_equilibrium_fixed_phase_volume(num_regression):
    '''
    Test a problem with H2O, NaCl, CO2, CaCO3 and Calcite 
    with fixed values of Phase Volumes 
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
  
    state = equilibrate(problem)
      
    output = stateDict(state)
      
    num_regression.check(output, 
                         default_tolerance=dict(atol=1e-7, rtol=1e-18))
  
  
 
#DO NOT TRY TO PUT ALL EQUILIBRIUM TEST IN A SINGUE TEST
#if one one of then fail the test will stop and you won't test all
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
def test_equilibrate_with_problem_setup(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives :
    - Equilibrium Problem 
    - Equilibrium Inverse Problem   
    '''   
    
    equilibriumState = ChemicalState();
    equilibriumState = equilibrate(problemSetup)
    result = stateDict(equilibriumState)
    num_regression.check(result)
  
 
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
def test_equilibritrate_with_problem_and_options(
        problemSetup,
        num_regression):
    '''
    checks equilibrate that receives:
    - Equilibrium Problem, Options
    - Equilibrium Inverse Problem, Options
    '''
    options = EquilibriumOptions()
    equilibriumState = equilibrate(problemSetup, options)
     
    output = stateDict(equilibriumState)
    num_regression.check(output)
     
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
def test_equilibrite_with_state(
        problemSetup,
        num_regression):
    '''
    checks equilibrate that receives:
    - state from Equilibrium Problem 
    - state from Equilibrium Inverse Problem
    '''
    #find a state based on the equilibrium problemSetup
    equilibriumState = equilibrate(problemSetup)
    #change pressure and temperature
    equilibriumState.setTemperature(problemSetup.temperature()+20)
    equilibriumState.setPressure(problemSetup.pressure()+20)
    #compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(equilibriumState)
     
    output = stateDict(equilibriumState)
    num_regression.check(output)
 
 
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
def test_equilibrate_with_state_partition(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives:
    - state and partition from Equilibrium Problem 
    - state and partition from Equilibrium Inverse Problem
    '''
    #find a state based on the equilibrium problemSetup  
    state = equilibrate(problemSetup)
    #change pressure and temperature
    state.setTemperature(problemSetup.temperature()+10)
    state.setPressure(problemSetup.pressure()+20)
     
    partition = Partition(problemSetup.system())
    #compute equilibrium for state with new temperature and pressure 
    equilibriumResult = equilibrate(state, partition)
       
    output = stateDict(state)
      
    num_regression.check(output)       
 
 
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
def test_equilibrate_with_state_option(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives:
    - state and options from Equilibrium Problem 
    - state and options from Equilibrium Inverse Problem
    '''
    #find a state based on the equilibrium problemSetup  
    state = equilibrate(problemSetup)
    #change pressure and temperature
    state.setTemperature(problemSetup.temperature()+10)
    state.setPressure(problemSetup.pressure()+20)
     
    options = EquilibriumOptions()
    #compute equilibrium for state with new temperature and pressure 
 
    equilibriumResult = equilibrate(state, options)
     
    output = stateDict(state)
      
    num_regression.check(output)       
 
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
def test_equilibrate_with_state_partition_option(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives:
    - state, partition and options from Equilibrium Problem 
    - state, partition and options from Equilibrium Inverse Problem
    '''
    #find a state based on the equilibrium problemSetup  
    state = equilibrate(problemSetup)
    #change pressure and temperature
    state.setTemperature(problemSetup.temperature()+10)
    state.setPressure(problemSetup.pressure()+20)
     
    options = EquilibriumOptions()
    partition = Partition(problemSetup.system())
     
    #compute equilibrium for state with new temperature and pressure 
    equilibriumResult = equilibrate(state, partition, options)
     
    output = stateDict(state)
      
    num_regression.check(output)       


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
def test_equilibrate_with_state_problemSetup(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives:
    - state and Equilibrium Problem 
    - state and Equilibrium Inverse Problem
    '''
    state = ChemicalState(problemSetup.system())
    
    equilibriumResult = equilibrate(state, problemSetup)    
    
    output = stateDict(state)
     
    num_regression.check(output) 
 
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
def test_equilibrate_with_state_problemSetup_option(
    problemSetup,
    num_regression):
    '''
    checks equilibrate that receives:
    - state, Equilibrium Problem and options 
    - state, Equilibrium Inverse Problem and options
    '''
    state = ChemicalState(problemSetup.system())
     
    options = EquilibriumOptions()
 
    equilibriumResult = equilibrate(state, problemSetup, options)
     
    output = stateDict(state)
      
    num_regression.check(output)
