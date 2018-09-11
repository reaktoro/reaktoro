import pytest
import io
import math
import numpy as np

from PyReaktoro import *
from pytest_regressions.plugin import num_regression
from _pytest.fixtures import fixture


def test_solve_state_T_P_be(num_regression):
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

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state,
                problem.temperature(),
                problem.pressure(),
                problem.elementAmounts())

    num_regression.check(stateDict(state))


def test_solve_state_problem(num_regression):
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

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state, problem)

    num_regression.check(stateDict(state))

def test_solve_state(num_regression):
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

    state = equilibrate(problem)
    state.setTemperature(problem.temperature()+30)
    state.setPressure(problem.pressure()+40)

    
    solver = EquilibriumSolver(system)

    solver.solve(state)

    num_regression.check(stateDict(state))


def test_approximate_state_T_P_be(num_regression):
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

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state,
                        problem.temperature(),
                        problem.pressure(),
                        problem.elementAmounts())

    num_regression.check(stateDict(state))

def test_approximate_state_problem(num_regression):
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

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state, problem)

    num_regression.check(stateDict(state))

def test_approximate_state(num_regression):
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

    state = equilibrate(problem)
    state.setTemperature(problem.temperature() + 30)
    state.setPressure(problem.pressure() + 40)

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
