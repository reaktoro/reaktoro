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
#   •
# -----------------------------------------------------------------------------

import reaktoro as rkt
import numpy as np
import math

# Define Thermofun database
db = rkt.ThermoFunDatabase ("psinagra-12-07")

# Define aqueous phase
solution = rkt.AqueousPhase(rkt.speciate("Al C Ca Cl Fe H K Mg N Na O P S Si U"))
solution.setActivityModel(rkt.chain(
    rkt.ActivityModelHKF(),
    rkt.ActivityModelDrummond("CO2")
))

# Define chemical system by providing database, aqueous phase, and minerals
system = rkt.ChemicalSystem(db, solution)

# Specify conditions to be satisfied at chemical equilibrium
specs = rkt.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

# Define equilibrium solver
solver = rkt.EquilibriumSolver(specs)

# Define temperature and pressure
T = 25.0 # in Celsius
P = 1.0 # in bar
#P = 0.01 # in bar

# Define conditions to be satisfied at chemical equilibrium
conditions = rkt.EquilibriumConditions(specs)
conditions.temperature(T, "celsius")
conditions.pressure(P, "bar")

#pHs = np.array([5, 6, 7, 8, 9, 10])
pHs = np.linspace(5, 10, num=41)
species_list = rkt.SpeciesList("UO2+2 UO2OH+ UO2(OH)2@ UO2CO3@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 "
                               "UO2(CO3)3-4 UO2(CO3)2-2 (UO2)2(OH)2+2 (UO2)3(OH)5+ (UO2)4(OH)7+")
bU = 0.0
nU = np.zeros((len(pHs), species_list.size()))

# Define initial equilibrium state with the following recipe:
# 1000 g H2O
# H3PO4@ 1e-5 mol
# CO2@ 0.1 g
 # UO2(SO4)@ 1e-5 mol.
state = rkt.ChemicalState(system)
state.set("H2O@"     , 1e3,  "g")
state.set("CO2@"     , 1e-1, "g")
state.set("H3PO4@"   , 1e-5, "mol")
state.set("UO2(SO4)@", 1e-5, "mol")

# Aqueous properties of the chemical state
aprops = rkt.AqueousProps(system)
props = rkt.ChemicalProps(system)

# Output the header of the table
# print("   pH  success" \
#       "          %n(UO2+2)" \
#       "         %n(UO2OH+)" \
#       "      %n(UO2(OH)2@)" \
#       "        %n(UO2CO3@)" \
#       "%n((UO2)3CO3(OH)3+)" \
#       "   %n((UO2)2(OH)+3)" \
#       "    %n(UO2(CO3)3-4)" \
#       "    %n(UO2(CO3)2-2)" \
#       "  %n((UO2)2(OH)2+2)" \
#       "   %n((UO2)3(OH)5+)" \
#       "   %n((UO2)4(OH)7+)")

for i in range(0, len(pHs)):

    # Set the value of pH for the current equilibrium calculations
    conditions.pH(pHs[i])

    # Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(state, conditions)

    if not res.optima.succeeded:
        nU[i, :] = math.nan * np.ones(species_list.size())
    else:
        # Update aqueous properties
        aprops.update(state)
        props.update(state)

        # Calculate U(VI) Speciation, %
        bU = props.elementAmount("U")[0]
        #print(f"{pHs[i]:6.2f}    {res.optima.succeeded}", end = " ")
        for j in range(0, species_list.size()):#species in species_list:
            nU[i, j] = state.speciesAmount(species_list[j].name())[0] / bU * 100
            #print(f"{nU[i, j]:18.2e}", end=' ')
        #print("")


import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9', 'C0', 'darkblue', 'darkgreen']

species_set_high = rkt.SpeciesList("UO2+2 UO2OH+ UO2CO3@ " \
                                   "UO2(CO3)3-4 UO2(CO3)2-2")
indices_high = [0, 1, 3, 6, 7]
species_set_low = rkt.SpeciesList("UO2(OH)2@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 " \
                                   "(UO2)2(OH)2+2 (UO2)3(OH)5+ (UO2)4(OH)7+")
indices_low = [2, 4, 5, 8, 9, 10]

plt.figure()
plt.xlabel("pH")
plt.ylabel("U(VI) Speciation, %")

for j in indices_high:
    plt.plot(pHs, nU[:, j], label=species_list[j].name(), color=colors[j])

plt.legend(loc="best")
plt.grid()
#plt.show()
plt.savefig(f'U-vs-pH-high-P-{P}.png', bbox_inches='tight')
plt.close()
plt.figure()
plt.xlabel("pH")
plt.ylabel("U(VI) Speciation, %")

for j in indices_low:
    plt.plot(pHs, nU[:, j], label=species_list[j].name(), color=colors[j])

plt.legend(loc="best")
plt.grid()
#plt.show()
plt.savefig(f'U-vs-pH-low-P-{P}.png', bbox_inches='tight')
plt.close()
