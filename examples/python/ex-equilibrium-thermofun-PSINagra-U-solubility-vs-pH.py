# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
#
# This library is free software you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

# -----------------------------------------------------------------------------
# 👏 Acknowledgements 👏
# -----------------------------------------------------------------------------
# This example was originally authored by:
#   • Svetlana Kyas (4 February 2022)
#
# and since revised by:
#   • G.D. Miron (1 April 2022)
# -----------------------------------------------------------------------------

from reaktoro import *
import numpy as np
import pandas as pd

# Define the Thermofun database
db = ThermoFunDatabase ("psinagra-12-07")

# Define the aqueous phase
solution = AqueousPhase(speciate("C Cl H O P S U"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define chemical system by providing database and aqueous phase
system = ChemicalSystem(db, solution)

# Specify conditions to be satisfied at chemical equilibrium
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()
specs.fugacity("CO2(g)")

# Define conditions to be satisfied at the chemical equilibrium state
conditions = EquilibriumConditions(specs)
conditions.temperature(25.0, "celsius")
conditions.pressure(1.0, "bar")
conditions.fugacity("CO2(g)", 0.01, "bar")

# Define the equilibrium solver and its options
solver = EquilibriumSolver(specs)
opts = EquilibriumOptions()

# Define initial equilibrium state
state = ChemicalState(system)
state.set("H2O@"     , 1e3,  "g")
state.set("UO2(OH)2@", 1e-5, "mol")

# Calculate the amount of uranium element
prop = ChemicalProps(state)
bU = prop.elementAmount("U")[0]

# Defined an auxiliary array of pH values
pHs = np.linspace(5, 10, num=61)

# Create list of species names, list of Species objects, and auxiliary amounts array
species_list_str = "UO2+2 UO2OH+ UO2(OH)2@ UO2CO3@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 " \
                   "UO2(CO3)3-4 UO2(CO3)2-2 (UO2)2(OH)2+2 (UO2)3(OH)5+ (UO2)4(OH)7+"
species_list = SpeciesList(species_list_str)
percentages = np.zeros(species_list.size())
amounts = np.zeros(species_list.size())

# Define dataframe to collect amount of the selected species
columns = ["pH"] + ["amount_" + name for name in species_list_str.split()] + ["perc_" + name for name in species_list_str.split()]
df = pd.DataFrame(columns=columns)

for pH in pHs:

    # Set the value of pH for the current equilibrium calculations
    conditions.pH(pH)

    # Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(state, conditions)

    # If the equilibrium calculations didn't succeed, continue to the next condition
    if not res.optima.succeeded: continue

    # Otherwise, calculate U(VI) Speciation, %
    for j in range(0, species_list.size()):#species in species_list:
        amounts[j] = float(state.speciesAmount(species_list[j].name()))
        percentages[j] = float(state.speciesAmount(species_list[j].name())) / bU * 100
    # Update dataframe with obtained values
    df.loc[len(df)] = np.concatenate([[pH], amounts, percentages])


import matplotlib.pyplot as plt
colors = ['teal', 'darkred', 'indigo', 'coral', 'rosybrown', 'steelblue', 'seagreen', 'palevioletred', 'darkred', 'darkkhaki', 'cadetblue', 'indianred']

plt.figure()
plt.xlabel("pH")
plt.ylabel("U(VI) Speciation, %")

species_list = species_list_str.split()
ax = df.plot(x="pH", y="perc_UO2+2", color=colors[0], label=species_list[0])

for species, color in zip(species_list[1:-1], colors[1:-1]):
    df.plot(x="pH", y="perc_"+species, ax=ax, color=color, label=species)
plt.legend(loc="best")
plt.grid()
plt.savefig(f'U-vs-pH.png', bbox_inches='tight')
plt.close()
