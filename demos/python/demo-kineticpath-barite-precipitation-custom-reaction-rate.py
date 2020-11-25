# Import reaktoro package
from reaktoro import *
import numpy as np
from math import *

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

# Initialize chemical system
system = ChemicalSystem(editor)

# Ionic strength function
I = ChemicalProperty.ionicStrength(system)

# Define barite mineral reaction and its parameters
eq_str_barite = "Barite = SO4-- + Ba++"
min_reaction_barite = editor.addMineralReaction("Barite") \
    .setEquation(eq_str_barite) \
    .addMechanism("logk = -8.6615 mol/(m2*s); Ea = 22 kJ/mol") \
    .setSpecificSurfaceArea(0.006, "m2/g")
# Create reaction for barite
reaction_barite = createReaction(min_reaction_barite, system)
reaction_barite.setName("Barite kinetics")
# Barite kinetic rate function
def rate_func_barite_shell(properties):


    # The number of chemical species in the system
    num_species = system.numSpecies()

    # Final mineral reaction rate using specified surface area
    res = ChemicalScalar(num_species)
    # Auxiliary variable for calculating mineral reaction rate
    f = ChemicalScalar(num_species, 1.0)

    # The universal gas constant (in units of kJ/(mol*K))
    R = 8.3144621e-3

    # The temperature and pressure of the system
    T = properties.temperature().val

    # Calculate the saturation index of the mineral
    lnK = reaction_barite.lnEquilibriumConstant(properties)
    lnQ = reaction_barite.lnReactionQuotient(properties)
    lnOmega = lnQ.val - lnK.val

    # Calculate the saturation index
    Omega = exp(lnOmega)

    # The composition of the chemical system
    n = properties.composition()

    # The initial composition of the chemical system
    n0 = reaction_barite.initialAmounts()

    # The index of the mineral
    imineral = system.indexSpeciesWithError(min_reaction_barite.mineral())

    # The number of moles of the mineral
    nm = n[imineral].val
    nm0 = n0[imineral]

    # Prevent negative mole numbers here for the solution of the ODEs
    nm = max(nm, 0.0)
    nm0 = max(nm0, 0.0)

    # Get H+ activity and ionic strength
    lna = properties.lnActivities().val
    i_h = system.indexSpeciesWithError("H+")

    activity_h = exp(lna[i_h])
    ionic_strength = I(properties).val

    # The molar mass of the mineral (in units of kg/mol)
    molar_mass = system.species(imineral).molarMass()

    # If (m <= 0) and (SRmin < 1) Then GoTo 240
    # If (SRmin = 1) Then GoTo 240
    if ((nm <= 0 and Omega < 1) or (Omega == 1)): # the is no way to precipitate further
        res = 0
    ssa = 0.006

    # If (SRmin > 1) Then GoTo 130
    if Omega > 1: # precipitation kinetics
        # If (m <= 1e-5) then GoTo 170
        if nm > 1e-5:
            kappa_pre = -2.18e-9 * exp(- 22 / R * (1.0/T - 1.0/298.15)) \
                        * pow(10, 0.6 * pow(ionic_strength, 0.5))
            res += f * ssa * nm * molar_mass * kappa_pre * (Omega - 1)
        else:
            # Set nucleation rate
            res = -1e-12 * f

    else: # dissolution kinetics
        kappa = 2.75e-8 * exp(- 25 / R * (1.0/T - 1.0/298.15)) \
                * pow(activity_h, 0.03) \
                * pow(10, 0.6 * pow(ionic_strength, 0.5))
        res += f * ssa * nm * molar_mass * pow(nm / nm0, 2 / 3) * kappa * (1 - pow(Omega, 0.2))

    return res
# Set custom barite kinetic rate
reaction_barite.setRate(rate_func_barite_shell)

# Initialize reaction system
reactions = ReactionSystem(system, [reaction_barite])

# Define the partition of the kinetic and equilibrium species
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

# Set initial amount of species in = the barite reaction
reaction_barite.setInitialAmounts(state_fw.speciesAmounts())

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
plt.savefig("barite-amount-with-customized-kinetic-rate.png")