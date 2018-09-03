import pytest
import io
from PyReaktoro import *
from pytest_regressions.plugin import file_regression

#TODO: try to use num_regression
def test_integration_demo_equilibrium_co2_brine(file_regression):
    database = Database(b"supcrt98.xml")

    editor = ChemicalEditor(database)#
    editor.addAqueousPhase(b"H2O NaCl CO2")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 0.1, b"mol")
    
    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open('state.txt', 'r') as file:
        result = file.read()
        file_regression.check(result)

#TODO: try to use num_regression        
def test_integration_demo_equilibrium_co2_brine_custom(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase([b"H2O(l)", b"H+", b"OH-", b"Na+", b"Cl-", b"HCO3-", b"CO2(aq)", b"CO3--"]) \
        .setActivityModelDrummondCO2()
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"]). \
        setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase(b"Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.setTemperature(60, b"celsius")
    problem.setPressure(300, b"bar")
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"CO2", 100, b"g")
    problem.add(b"NaCl", 1, b"mol")

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open('state.txt', 'r') as file:
        result = file.read()
        file_regression.check(result)

#TODO: try to use num_regression        
def test_integration_demo_equilibrium_co2_brine_using_equilibriumsolver(file_regression):
   # Create an object of Database class to use in the initialization of the
   # chemical system.
   database = Database("supcrt98.xml")

   # Define which phases and species the chemical system should have using a
   # ChemicalEditor object.
   editor = ChemicalEditor(database)
   editor.addAqueousPhase(b"H2O NaCl CO2")
   editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
   editor.addMineralPhase(b"Halite")

   # Initialize the chemical system using the definition provided to the
   # chemical editor.
   system = ChemicalSystem(editor)

   # Define an equilibrium problem with given temperature, pressure, and
   # amounts of compounds.
   problem = EquilibriumProblem(system)
   problem.setTemperature(60, b"celsius")
   problem.setPressure(300, b"bar")
   problem.add(b"H2O", 1, b"kg")
   problem.add(b"CO2", 100, b"g")
   problem.add(b"NaCl", 0.1, b"mol")

   # Get the temperature, pressure, and mole amounts of the elements.
   T = problem.temperature()
   P = problem.pressure()
   b = problem.elementAmounts()
    
   # Create an object of EquilibriumSolver class that can be reused many times.
   solver = EquilibriumSolver(system)
   
   # Create an object of ChemicalState class to store the equilibrium
   # state of the system.
   state = ChemicalState(system)
   
   # Solve the equilibrium state with given (T, P, b) inputs.
   solver.solve(state, T, P, b)
    
   # Print the calculated chemical equilibrium state.
   state.output("state.txt")
    
   file = io.open('state.txt', 'r')
   with io.open('state.txt', 'r') as file:
       result = file.read()
       file.close()
       
   # Calculate the new equilibrium state when temperature is increased.
   # Use the previous equilibrium state as an initial guess for improved
   # performance.
   solver.solve(state, T + 10.0, P, b)

   # Print the new calculated chemical equilibrium state.
   state.output("state.txt")
   with io.open('state.txt', 'r') as file:
       result = result + file.read()
       file.close()
   
   file_regression.check(result)

#TODO: try to use num_regression
@pytest.mark.xfail(reason='RES-9')
def test_demo_equilibrium_fixed_alkalinity(file_regression):
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

    state.output("state.txt")
    
    with io.open("state.txt", 'r') as file:
        result = file.read()
        file_regression.check(result)
