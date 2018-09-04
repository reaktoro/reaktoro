import pytest
import io
from numpy import *
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
@pytest.mark.xfail(reason='RES-9')       
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
    
    state.output("state2.txt")
    
    with io.open('state2.txt', 'r') as file:
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

def test_demo_equilibrium_fixed_amount(file_regression):
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

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)
        
def test_demo_equilibrium_fixed_ph_activity_amount(file_regression):
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

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)    

def test_demo_equilibrium_fixed_ph_co2(file_regression):
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

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)

def test_demo_equilibrium_fixed_ph_hcl_naoh(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O Na Cl Ca C")
    editor.addMineralPhase(b"Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 0.1, b"mol")
    problem.pH(8.0, b"HCl", b"NaOH")
    problem.fixSpeciesAmount(b"Calcite", 1, b"mol")

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)

@pytest.mark.xfail(reason='RES-10')
def test_demo_equilibrium_fixed_phase_volume(file_regression):
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
    
    state.output("state.txt")

    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)    

def test_demo_equilibrium_iron_nh3(file_regression):
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhase(b"H2O Fe(OH)2 Fe(OH)3 NH3")
    editor.addGaseousPhase(b"NH3(g)")
    editor.addMineralPhase(b"Magnetite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"Fe(OH)2", 1, b"mol")
    problem.add(b"Fe(OH)3", 2, b"mol")
    problem.add(b"NH3", 1, b"mmol")

    state = equilibrate(problem)
    
    state.output("state.txt")
    
    with io.open("state.txt", "r") as file:
        result = file.read()
        file_regression.check(result)
        
    
@pytest.mark.xfail(reason='RES-10')
def test_demo_equilibriumpath_calcite_hcl(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O Ca C Cl")
    editor.addMineralPhase(b"Calcite")

    system = ChemicalSystem(editor)

    problem1 = EquilibriumProblem(system)
    problem1.setTemperature(30.0, b"celsius")
    problem1.setPressure(1.0, b"bar")
    problem1.add(b"H2O", 1, b"kg")
    problem1.add(b"CaCO3", 1, b"g")

    problem2 = EquilibriumProblem(system)
    problem2.setTemperature(30.0, b"celsius")
    problem2.setPressure(1.0, b"bar")
    problem2.add(b"H2O", 1, b"kg")
    problem2.add(b"CaCO3", 1, b"g")
    problem2.add(b"HCl", 1, b"mmol")

    state1 = equilibrate(problem1)
    state2 = equilibrate(problem2)

    path = EquilibriumPath(system)

    plot1 = path.plot()
    plot1.x(b"elementAmount(Cl units=mmol)")
    plot1.y(b"pH")
    plot1.xlabel(b"HCl [mmol]")
    plot1.ylabel(b"pH")
    plot1.showlegend(False)

    plot2 = path.plot()
    plot2.x(b"elementAmount(Cl units=mmol)")
    plot2.y(b"elementMolality(Ca units=mmolal)", b"Ca")
    plot2.xlabel(b"HCl [mmol]")
    plot2.ylabel(b"Concentration [mmolal]")
    plot2.legend(b"right center")
    
    plot3 = path.plot()
    plot3.x(b"elementAmount(Cl units=mmol)")
    plot3.y(b"speciesMolality(CO2(aq) units=mmolal)", b"CO2(aq)")
    plot3.y(b"speciesMolality(CO3-- units=mmolal)", b"CO3--")
    plot3.xlabel(b"HCl [mmol]")
    plot3.ylabel(b"Concentration [mmolal]")
    plot3.legend(b"right top")

    plot4 = path.plot()
    plot4.x(b"elementAmount(Cl units=mmol)")
    plot4.y(b"speciesMass(Calcite units=g)", b"Calcite")
    plot4.xlabel(b"HCl [mmol]")
    plot4.ylabel(b"Mass [g]")
    
    output = path.output()
    output.filename(b"result.txt")
    output.add(b"elementAmount(Cl units=mmol)", b"Cl [mmol]")
    output.add(b"elementMolality(Ca units=mmolal)", b"Ca [mmolal]")
    output.add(b"pH")
    output.add(b"speciesMass(Calcite units=g)")

    path.solve(state1, state2)
    
    with io.open('plot0.dat', 'r') as file:
       result = file.read()
       file.close()
    
    with io.open('plot1.dat', 'r') as file:
       result = result + file.read()
       file.close()
    
    with io.open('plot2.dat', 'r') as file:
       result = result + file.read()
       file.close()
       
    with io.open('plot3.dat', 'r') as file:
       result = result + file.read()
       file.close()
    
    with io.open('result.txt', 'r') as file:
       result = result + file.read()
       file.close()
       
    file_regression.check(result)


def test_demo_equilibriumpath_co2(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H O C Na Cl")

    system = ChemicalSystem(editor)

    problem1 = EquilibriumProblem(system)
    problem1.add(b"H2O", 1, b"kg")
    problem1.add(b"CO2", 0.5, b"mol")
    problem1.add(b"HCl", 1, b"mol")

    problem2 = EquilibriumProblem(system)
    problem2.add(b"H2O", 1, b"kg")
    problem2.add(b"CO2", 0.5, b"mol")
    problem2.add(b"NaOH", 2, b"mol")

    state1 = equilibrate(problem1)
    state2 = equilibrate(problem2)

    path = EquilibriumPath(system)

    plot = path.plot()
    plot.x(b"pH")
    plot.y(b"speciesMolality(HCO3-)", b"HCO@_3^-")
    plot.y(b"speciesMolality(CO2(aq))", b"CO_2(aq)")
    plot.y(b"speciesMolality(CO3--)", b"CO@_3^{2-")
    plot.xlabel(b"pH")
    plot.ylabel(b"Concentration [molal]")
    plot.yformat(b"%g")
    plot.legend(b"left center Left reverse")
    plot.name(b"CO2Pathplot")

    output = path.output()
    output.filename(b"result.txt")
    output.add(b"t")
    output.add(b"pH")
    output.add(b"speciesMolality(HCO3-)", b"HCO3- [molal]")
    output.add(b"speciesMolality(CO2(aq))", b"CO2(aq) [molal]")
    output.add(b"speciesMolality(CO3--)", b"CO3-- [molal]")

    path.solve(state1, state2)
    
        
    with io.open('CO2Pathplot.dat', 'r') as file:
       result = file.read()
       file.close()
        
    with io.open('result.txt', 'r') as file:
       result = result + file.read()
       file.close()
    
    file_regression.check(result)
    
def test_demo_kineticpath_calcite_hcl(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H2O HCl CaCO3")
    editor.addMineralPhase(b"Calcite")

    editor.addMineralReaction(b"Calcite") \
        .setEquation(b"Calcite = Ca++ + CO3--") \
        .addMechanism(b"logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")

    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)

    partition = Partition(system)
    partition.setKineticPhases([b"Calcite"])

    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"HCl", 1, b"mmol")

    state0 = equilibrate(problem)

    state0.setSpeciesMass(b"Calcite", 100, b"g")

    path = KineticPath(reactions)
    path.setPartition(partition)
    
    plot1 = path.plot()
    plot1.x(b"time(units=minute)")
    plot1.y(b"elementMolality(Ca units=mmolal)", b"Ca")
    plot1.xlabel(b"Time [minute]")
    plot1.ylabel(b"Concentration [mmolal]")
    plot1.legend(b"right center")
    plot1.name(b"a")
    
    plot2 = path.plot()
    plot2.x(b"time(units=minute)")
    plot2.y(b"phaseMass(Calcite units=g)", b"Calcite")
    plot2.xlabel(b"Time [minute]")
    plot2.ylabel(b"Mass [g]")
    plot2.name(b"b")

    path.solve(state0, 0, 5, b"minute")
    
            
    with io.open('a.dat', 'r') as file:
       result = file.read()
       file.close()
        
    with io.open('b.dat', 'r') as file:
       result = result + file.read()
       file.close()
    
    file_regression.check(result)
    

def test_demo_kineticpath_carbonates_co2(file_regression):
    editor = ChemicalEditor()
    editor.addAqueousPhase(b"H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase([b"H2O(g)", b"CO2(g)"])
    editor.addMineralPhase(b"Calcite")
    editor.addMineralPhase(b"Magnesite")
    editor.addMineralPhase(b"Dolomite")
    editor.addMineralPhase(b"Halite")

    editor.addMineralReaction(b"Calcite") \
        .setEquation(b"Calcite = Ca++ + CO3--") \
        .addMechanism(b"logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")

    editor.addMineralReaction(b"Magnesite") \
        .setEquation(b"Magnesite = Mg++ + CO3--") \
        .addMechanism(b"logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism(b"logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, b"cm2/g")

    editor.addMineralReaction(b"Dolomite") \
        .setEquation(b"Dolomite = Ca++ + Mg++ + 2*CO3--") \
        .addMechanism(b"logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
        .addMechanism(b"logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
        .setSpecificSurfaceArea(10, b"cm2/g")

    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)

    partition = Partition(system)
    partition.setKineticSpecies([b"Calcite", b"Magnesite", b"Dolomite"])

    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add(b"H2O", 1, b"kg")
    problem.add(b"NaCl", 1, b"mol")
    problem.add(b"CO2", 1, b"mol")

    state0 = equilibrate(problem)

    state0.setSpeciesMass(b"Calcite", 100, b"g")
    state0.setSpeciesMass(b"Dolomite", 50, b"g")

    path = KineticPath(reactions)
    path.setPartition(partition)

    plot0 = path.plot()
    plot0.x(b"time(units=hour)")
    plot0.y(b"pH")
    plot0.xlabel(b"Time [hour]")
    plot0.ylabel(b"pH")
    plot0.showlegend(False)

    plot1 = path.plot()
    plot1.x(b"time(units=hour)")
    plot1.y(b"elementMolality(Ca)", b"Ca")
    plot1.y(b"elementMolality(Mg)", b"Mg")
    plot1.xlabel(b"Time [hour]")
    plot1.ylabel(b"Concentration [molal]")
    plot1.legend(b"right center")
    plot1.name("c")

    plot2 = path.plot()
    plot2.x(b"time(units=hour)")
    plot2.y(b"phaseMass(Calcite units=grams)", b"Calcite")
    plot2.xlabel(b"Time [hour]")
    plot2.ylabel(b"Mass [g]")
    plot2.name("d")
    
    plot3 = path.plot()
    plot3.x(b"time(units=hour)")
    plot3.y(b"phaseMass(Dolomite units=grams)", b"Dolomite")
    plot3.xlabel(b"Time [hour]")
    plot3.ylabel(b"Mass [g]")
    plot3.name("e")
    
    path.solve(state0, 0, 25, b"hours")
    
    state0.output("state3.txt")

    with io.open('c.dat', 'r') as file:
       result = file.read()
       file.close()
        
    with io.open('d.dat', 'r') as file:
       result = result + file.read()
       file.close()
       
    with io.open('e.dat', 'r') as file:
       result = result + file.read()
       file.close()
    
    with io.open('state3.txt', 'r') as file:
       result = result + file.read()
       file.close()
       
    file_regression.check(result) 

    
def test_demo_thermodynamic_properties_of_species_using_supcrt_database(file_regression):
    # Create a Database object loaded with SUPCRT98 database file.
    database = Database("supcrt98.xml")

    # Create a Thermo object for thermodynamic property calculations
    thermo = Thermo(database)

    # Create shorter alias for the methods in Thermo that calculate
    # standard partial molar Gibbs energy and enthalpy of species
    evalG0 = thermo.standardPartialMolarGibbsEnergy
    evalH0 = thermo.standardPartialMolarEnthalpy

    # Create a list of species names, as found in the database
    species = ['H2O(l)', 'HCO3-', 'CO2(aq)', 'CO2(g)', 'Calcite']

    # Create a numeric array of temperature values (in units of K)
    temperatures = array([25, 50, 75, 100, 200, 300]) + 273.15

    # Create pressure variable (in units of Pa)
    P = 100.e+5

    # Create Python dictionaries containing the standard partial
    # molar Gibbs energy and enthalpy for each species
    G0 = {}
    H0 = {}

    # Calculate the standard chemical potentials of the species at
    # 100 bar and at each temperature in the array of temperatures
    for name in species:
        G0[name] = array([evalG0(T, P, name).val for T in temperatures])
        H0[name] = array([evalH0(T, P, name).val for T in temperatures])

        # For each species, create a file containing three columns:
        # 1st column: temperature (in units of K)
        # 2nd column: standard partial molar gibbs energy (in units of J/mol)
        # 3rd column: standard partial molar enthalpy (in units of J/mol)
    with io.open('calculated-standard-species-properties.txt', 'w') as file:
        for name in species:
            file.write(repr(name) + "\n") 
            file.write("T(K), G0(J/mol), H0(J/mol)\n")
            for Tval, G0val, H0val in zip(temperatures, G0[name], H0[name]):
                file.write("{0}, {1}, {2}\n".format(Tval, G0val, H0val))
    
    with io.open('calculated-standard-species-properties.txt', 'r') as file:
        result = file.read()
        file.close()
    
    file_regression.check(result)  

def test_demo_thermo_properties(file_regression):
    # Create a Database object loaded with SUPCRT database file.
    database = Database("supcrt98.xml")

    # Create a Thermo object for thermodynamic property calculations
    thermo = Thermo(database)

    # Define two reactions for which we'll calculate their log(K)
    reaction1 = 'CO2(g) + H2O(l) = HCO3- + H+'
    reaction2 = 'Calcite + H+ = Ca++ + HCO3-'

    # Define pressure as 100 bar for the calculations (converted to Pa)
    P = 100.0e5

    # Create an array with temperature values from 25 to 300 C (converted to K)
    x = linspace(25.0, 300.0, 50) + 273.15

    # Calculate the log(K) of the reactions at those temperature points and 100 bar
    y1 = [thermo.logEquilibriumConstant(T, P, reaction1).val for T in x]
    y2 = [thermo.logEquilibriumConstant(T, P, reaction2).val for T in x]

    x = x - 273.15

    with io.open('calculated-properties.txt', 'w') as file:
        file.write(reaction1 + "\n") 
        file.write("T, P\n")
        for i in range(0,len(x)):
            file.write('{0} {1} {2}\n'.format(x[i], y1[i], y2[i]))
    
    with io.open('calculated-properties.txt', 'r') as file:
        result = file.read()
        file.close()
    
    file_regression.check(result)    
     
#TODO - add a test based on demo-water-properties.py
def test_demo_water_properties(file_regression):
    T, P = 298.15, 1e5

    wts = waterThermoStateWagnerPruss(T, P, StateOfMatter.Liquid)
    
    with io.open('water-properties.txt', 'w') as file:
        file.write("temperature.val = {0}".format(wts.temperature.val))
        file.write("temperature.ddT = {0}".format(wts.temperature.ddT))
        file.write("temperature.ddP = {0}".format(wts.temperature.ddP))
        file.write("volume.val = {0}".format(wts.volume.val))
        file.write("volume.ddT = {0}".format(wts.volume.ddT))
        file.write("volume.ddP = {0}".format(wts.volume.ddP))
        file.write("entropy.val = {0}".format(wts.entropy.val))
        file.write("entropy.ddT = {0}".format(wts.entropy.ddT))
        file.write("entropy.ddP = {0}".format(wts.entropy.ddP))
        file.write("helmholtz.val = {0}".format(wts.helmholtz.val))
        file.write("helmholtz.ddT = {0}".format(wts.helmholtz.ddT))
        file.write("helmholtz.ddP = {0}".format(wts.helmholtz.ddP))
        file.write("internal_energy.val = {0}".format(wts.internal_energy.val))
        file.write("internal_energy.ddT = {0}".format(wts.internal_energy.ddT))
        file.write("internal_energy.ddP = {0}".format(wts.internal_energy.ddP))
        file.write("enthalpy.val = {0}".format(wts.enthalpy.val))
        file.write("enthalpy.ddT = {0}".format(wts.enthalpy.ddT))
        file.write("enthalpy.ddP = {0}".format(wts.enthalpy.ddP))
        file.write("gibbs.val = {0}".format(wts.gibbs.val))
        file.write("gibbs.ddT = {0}".format(wts.gibbs.ddT))
        file.write("gibbs.ddP = {0}".format(wts.gibbs.ddP))
        file.write("cv.val = {0}".format(wts.cv.val))
        file.write("cv.ddT = {0}".format(wts.cv.ddT))
        file.write("cv.ddP = {0}".format(wts.cv.ddP))
        file.write("cp.val = {0}".format(wts.cp.val))
        file.write("cp.ddT = {0}".format(wts.cp.ddT))
        file.write("cp.ddP = {0}".format(wts.cp.ddP))
        file.write("density.val = {0}".format(wts.density.val))
        file.write("density.ddT = {0}".format(wts.density.ddT))
        file.write("density.ddP = {0}".format(wts.density.ddP))
        file.write("densityT.val = {0}".format(wts.densityT.val))
        file.write("densityT.ddT = {0}".format(wts.densityT.ddT))
        file.write("densityT.ddP = {0}".format(wts.densityT.ddP))
        file.write("densityP.val = {0}".format(wts.densityP.val))
        file.write("densityP.ddT = {0}".format(wts.densityP.ddT))
        file.write("densityP.ddP = {0}".format(wts.densityP.ddP))
        file.write("densityTT.val = {0}".format(wts.densityTT.val))
        file.write("densityTT.ddT = {0}".format(wts.densityTT.ddT))
        file.write("densityTT.ddP = {0}".format(wts.densityTT.ddP))
        file.write("densityTP.val = {0}".format(wts.densityTP.val))
        file.write("densityTP.ddT = {0}".format(wts.densityTP.ddT))
        file.write("densityTP.ddP = {0}".format(wts.densityTP.ddP))
        file.write("densityPP.val = {0}".format(wts.densityPP.val))
        file.write("densityPP.ddT = {0}".format(wts.densityPP.ddT))
        file.write("densityPP.ddP = {0}".format(wts.densityPP.ddP))
        file.write("pressure.val = {0}".format(wts.pressure.val))
        file.write("pressure.ddT = {0}".format(wts.pressure.ddT))
        file.write("pressure.ddP = {0}".format(wts.pressure.ddP))
        file.write("pressureT.val = {0}".format(wts.pressureT.val))
        file.write("pressureT.ddT = {0}".format(wts.pressureT.ddT))
        file.write("pressureT.ddP = {0}".format(wts.pressureT.ddP))
        file.write("pressureD.val = {0}".format(wts.pressureD.val))
        file.write("pressureD.ddT = {0}".format(wts.pressureD.ddT))
        file.write("pressureD.ddP = {0}".format(wts.pressureD.ddP))
        file.write("pressureTT.val = {0}".format(wts.pressureTT.val))
        file.write("pressureTT.ddT = {0}".format(wts.pressureTT.ddT))
        file.write("pressureTT.ddP = {0}".format(wts.pressureTT.ddP))
        file.write("pressureTD.val = {0}".format(wts.pressureTD.val))
        file.write("pressureTD.ddT = {0}".format(wts.pressureTD.ddT))
        file.write("pressureTD.ddP = {0}".format(wts.pressureTD.ddP))
        file.write("pressureDD.val = {0}".format(wts.pressureDD.val))
        file.write("pressureDD.ddT = {0}".format(wts.pressureDD.ddT))
        file.write("pressureDD.ddP = {0}".format(wts.pressureDD.ddP))
 
    with io.open('water-properties.txt', 'r') as file:
        result = file.read()
        file.close()
    
    file_regression.check(result)  