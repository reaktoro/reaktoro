# Import reaktoro package
from reaktoro import *
import numpy as np

# Define time related constants
second = 1
minute = 60 * second
hour = 60 * minute
day = 60 * hour

# Thermodynamic conditions
T = 60.0        # temperature (in units of Celsius)
P = 200.0       # pressure (in units of atm)
water_kg = 1.00 # water mass

# Define database
db = Database('supcrt07.xml')

# Fetch Debye-Huckel activity model parameters
dhModel = DebyeHuckelParams()
dhModel.setPHREEQC()

# Define different phases of the chemical system
editor = ChemicalEditor(db)
editor.addAqueousPhaseWithElements("H Cl S O Ba Ca Sr Na K Mg C Si"). \
    setChemicalModelDebyeHuckel(dhModel)
editor.addMineralPhase('Barite')    # BaSiO4

# Define barite mineral reaction and its parameters
eq_str_barite = "Barite = SO4-- + Ba++"
min_reaction_barite = editor.addMineralReaction("Barite") \
    .setEquation(eq_str_barite) \
    .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol") \
    .setSpecificSurfaceArea(0.006, "m2/g")

# Initialize chemical system
system = ChemicalSystem(editor)
# Initialize reaction system
reactions = ReactionSystem(editor)

partition = Partition(system)
partition.setKineticSpecies(["Barite"])

# Defining the formation water chemical state
problem_fw = EquilibriumInverseProblem(system)
problem_fw.setTemperature(T, "celsius")
problem_fw.setPressure(P, "atm")
problem_fw.add('H2O', water_kg, 'kg')
problem_fw.add("SO4", 10 * water_kg, "ug")
problem_fw.add("Ca", 995 * water_kg, "mg")
problem_fw.add("Ba", 995 * water_kg, "mg")
problem_fw.add("Sr", 105 * water_kg, "mg")
problem_fw.add("Na", 27250 * water_kg, "mg")
problem_fw.add("K", 1730 * water_kg, "mg")
problem_fw.add("Mg", 110 * water_kg, "mg")
problem_fw.add("Cl", 45150 * water_kg, "mg")
problem_fw.add("HCO3", 1980 * water_kg, "mg")
problem_fw.pH(7.0, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_fw = equilibrate(problem_fw)
state_fw.setSpeciesAmount("Barite", 0.1, "mcmol")
state_fw.scaleVolume(1.0, "m3")

# Defining the seawater chemical state
problem_sw = EquilibriumInverseProblem(system)
problem_sw.setTemperature(T, "celsius")
problem_sw.setPressure(P, "atm")
problem_sw.add('H2O', water_kg, 'kg')
problem_sw.add("SO4--", 2710 * water_kg, "mg")
problem_sw.add("Ca++", 411 * water_kg, "mg")
problem_sw.add("Ba++", 0.01 * water_kg, "mg")
problem_sw.add("Sr++", 8 * water_kg, "mg")
problem_sw.add("Na+", 10760 * water_kg, "mg")
problem_sw.add("K+", 399 * water_kg, "mg")
problem_sw.add("Mg++", 1290 * water_kg, "mg")
problem_sw.add("Cl-", 19350 * water_kg, "mg")
problem_sw.add("HCO3-", 142 * water_kg, "mg")
problem_sw.pH(8.1, "HCl", "NaOH")

# Calculate the equilibrium states for the initial conditions
state_sw = equilibrate(problem_sw)
state_sw.scaleVolume(1.0, "m3")

# Solving kinetic barite precipitation problem
t0, tfinal = 0.0, 1.0
result_file_name = "kineticpath-barite-precipitation-tfinal-" + str(tfinal*day) + ".txt"
path = KineticPath(reactions, partition)

# Define the characteristics of the chemical state to be output
output = path.output()
output.filename(result_file_name)
output.add("time(units=s)")
output.add("speciesAmount(Barite units=mol)", "Barite")
output.filename(result_file_name)

# Perturb the formation water state with seawater state
state_fw = state_fw + state_sw
# Solve the kinetic path problem
path.solve(state_fw, t0, tfinal, "days")

# Plotting of the results of kinetic path calculation
filearray = np.loadtxt(result_file_name, skiprows=1) # load data from the file skipping the one row
data = filearray.T  # transpose the matrix with data
[time_indx, barite_indx] = np.arange(0, 2)

# To visually analyze the obtained reaction path, we export
# [bokeh](https://docs.bokeh.org/en/latest/docs/gallery.html#standalone-examples) python plotting package.

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.set_loglevel("critical")

plt.figure()
plt.plot(data[time_indx], data[barite_indx], label="Barite", color="C1")
plt.legend(loc='center right')
plt.xlabel("Time [s]")
plt.ylabel("Amount [mol]")
plt.savefig("barite-amount.png")