import pytest
import io
import math
import numpy as np

from PyReaktoro import *
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture

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

@pytest.fixture
def problemSetupH2O_CO2_NaCl_Halite_60C_300P():
    '''
    Build a problem with 1kg of H2O, 100g of CO2 and 0.1mol of NaCl 
    at 60ºC and 300 bar 
    '''
    database = Database(b"supcrt98.xml")
    
    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O NaCl CO2")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Halite")
    
    system = ChemicalSystem(editor)
    
    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    
    return problem

@pytest.fixture
def problemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P():
    '''
    Build a problem with H2O, H+, Na+, Cl-, HCO3-, CO2(aq), CO3-- and 
    Halite at 60ºC and 300 bar, defining their chemical species 
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase([b"H2O(l)", b"H+", b"OH-", b"Na+", b"Cl-", b"HCO3-", b"CO2(aq)", b"CO3--"]) \
        .setActivityModelDrummondCO2()
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"]). \
        setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase(b"Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 1, b"mol")
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    
    return problem
    
@pytest.fixture
def problemSetupH2O_NaCl_CaCO3_CalcilteFixedMass():
    '''
    Build a problem with H2O, NaCL, CaCO3, Calcite with fixed
    species mass and amount  
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
    
    return problem

@pytest.fixture
def problemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with fixed
    species amount, activity and defined pH  
    '''
    
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O Na Cl Ca Mg C")
    editor.addGaseousPhase(b"H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.add(b"CaCl2", 2, b"mmol")
    problem.add(b"MgCl2", 4, b"mmol")
    problem.pH(3.0, b"HCl")
    problem.fixSpeciesAmount(b"CO2(g)", 1.0, b"mol")
    problem.fixSpeciesActivity(b"O2(g)", 0.20)
    
    return problem

@pytest.fixture
def problemSetupH_O_Na_Cl_Ca_Mg_CpH():
    '''
    Build a problem with H, Na, Cl, Ca, Mg, C with defined pH  
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O Na Cl Ca Mg C")
    editor.addGaseousPhase(b"H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.add(b"CaCl2", 2, b"mmol")
    problem.add(b"MgCl2", 4, b"mmol")
    problem.pH(4.0, b"CO2")
    
    return problem

@pytest.fixture
def problemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount():
    '''
    Build a problem with H, O, Na, Cl, Ca, C and Calcite with defined pH 
    and fixed species amount  
    '''
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O Na Cl Ca C")
    editor.addMineralPhase(b"Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.pH(8.0, b"HCl", b"NaOH")
    problem.fixSpeciesAmount(b"Calcite", 1, b"mol")

    return problem
    
    
@pytest.mark.parametrize('problemSetup',
    [
        (
            pytest.lazy_fixture('problemSetupH2O_CO2_NaCl_Halite_60C_300P')
        ),
        (
            pytest.lazy_fixture('problemSetupH2O_CO2_NaCl_Halite_dissolved_60C_300P')
        ),
        (
            pytest.lazy_fixture('problemSetupH2O_NaCl_CaCO3_CalcilteFixedMass')
        ),
        (
            pytest.lazy_fixture('problemSetupH_O_Na_Cl_Ca_Mg_CFixedAmountAndActivity')
        ),
        (
            pytest.lazy_fixture('problemSetupH_O_Na_Cl_Ca_Mg_CpH')
        ),
        (
            pytest.lazy_fixture('problemSetupH_O_Na_Cl_Ca_C_CalcitepHFixedAmount')
        )
    ],
    ids=['H2O, CO2, NaCl and Halite at 60C and 300 bar',
         'H2O, CO2, NaCl and Halite already dissolved at 60C and 300 bar',
         'H2O, CO2, NaCl, CaCO3, Calcite with fixed species mass and amounts',
         'H, O, Na, Cl, Ca, Mg, C with fixed amount, activity and pH',
         'H, O, Na, Cl, Ca, Mg, C with defined pH',
         'H, O, Na, Cl, Ca, C, Calcite with defined pH and fixed amount'
         ]
    )
def test_equilibrium_calculation(
    problemSetup,
    num_regression):
    '''
    checks equilibrium state calculation for given problem setup 
    '''    
    
    state = equilibrate(problemSetup)
    
    output = stateDict(state)
    
    num_regression.check(output, 
                         default_tolerance=dict(atol=1e-7, rtol=1e-18))



def test_equilibrium_calculation_using_equilibriumsolver(
    problemSetupH2O_CO2_NaCl_Halite_60C_300P,
    num_regression):
    '''
    calculates multiples equilibrium state for the given setup problem
    '''    
    system = problemSetupH2O_CO2_NaCl_Halite_60C_300P.system()
    
    solver = EquilibriumSolver(system)
    
    state = ChemicalState(system)
    
    solver.solve(state, 
                problemSetupH2O_CO2_NaCl_Halite_60C_300P.temperature(),
                problemSetupH2O_CO2_NaCl_Halite_60C_300P.pressure(),
                problemSetupH2O_CO2_NaCl_Halite_60C_300P.elementAmounts())
    
    outputState1 = stateDict(state)
    
    num_regression.check(outputState1, 
                        basename="test_equilibrium_calculation_H2O_NaCl_CO2_using_equilibriumsolver_state1.txt",
                        tolerances=None, 
                        default_tolerance=dict(atol=1e-7, rtol=1e-18))
     
    solver.solve(state,
                 problemSetupH2O_CO2_NaCl_Halite_60C_300P.temperature() + 10.0,
                 problemSetupH2O_CO2_NaCl_Halite_60C_300P.pressure(),
                 problemSetupH2O_CO2_NaCl_Halite_60C_300P.elementAmounts())
    
    outputState2 = stateDict(state)
    
    num_regression.check(outputState2, 
                        basename="test_equilibrium_calculation_H2O_NaCl_CO2_using_equilibriumsolver_state2.txt",
                        tolerances=None,  
                        default_tolerance=dict(atol=1e-7, rtol=1e-18))
    
@pytest.mark.xfail(reason='RES-9')
def test_demo_equilibrium_fixed_alkalinity(num_regression):
    '''
    Build a problem with H2O, NaCl, CO2, CaCO3 and Calcite 
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
    