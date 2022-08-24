# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

# -----------------------------------------------------------------------------
# 👏 Acknowledgements 👏
# -----------------------------------------------------------------------------
# This example was originally authored by:
#   • Svetlana Kyas (13 July 2022)
# -----------------------------------------------------------------------------

# https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-64.htm

from reaktoro import *
import numpy as np
import math

# Initialize a thermodynamic database
db = PhreeqcDatabase("phreeqc.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Cl Na K Mg S Si"))
aqueousphase.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Gypsum
# 	CaSO4:2H2O = Ca+2 + SO4-2 + 2 H2O
# 	-log_k	-4.58
# 	-delta_h -0.109 kcal
# 	-analytic	68.2401	0.0	-3221.51	-25.0627
# 	-analytical_expression  93.7  5.99E-03  -4e3  -35.019 # better fits the appendix data of Appelo, 2015, AG 55, 62
# 	-Vm 73.9 # 172.18 / 2.33  (Vm H2O = 13.9 cm3/mol)
# Anhydrite
# 	CaSO4 = Ca+2 + SO4-2
# 	-log_k	-4.36
# 	-delta_h -1.710 kcal
# 	-analytic  84.90  0  -3135.12  -31.79 # 50 - 160oC, 1 - 1e3 atm, anhydrite dissolution, Blount and Dickson, 1973, Am. Mineral. 58, 323.
# 	-Vm 46.1 # 136.14 / 2.95
# Create mineral phases
minerals = MineralPhases("Gypsum Anhydrite")

# Construct the chemical system
system = ChemicalSystem(db, aqueousphase, minerals)
props = ChemicalProps(system)
aprops = AqueousProps(system)

# Create the equilibrium solver
solver = EquilibriumSolver(system)

# Initial values of the carbonates and mass of water
n0Gypsum = 10.0
n0Anhydrite = 10.0
water_kg = 1.0

# Function to calculate equilibrium of the carbonates and seawater
def minerals_in_purewater(T, P):

    state = ChemicalState(system)
    state.temperature(T, "celsius")
    state.pressure(P, "atm")
    state.set("H2O", water_kg, "kg")
    # Add minerals
    state.set("Gypsum", n0Gypsum, "mol")
    state.set("Anhydrite", n0Anhydrite, "mol")

    # Calculate chemical state corresponding to the seawater
    res = solver.solve(state)

    # Equilibrate the seawater with carbonates
    props.update(state)
    aprops.update(state)

    # Throw exception if the equilibrium couldn't be found
    if not res.optima.succeeded:
        nGypsum = math.nan
        nAnhydrite = math.nan
        mCa = math.nan
    else:
        # Fetch values of the specified species
        nGypsum    = float(state.speciesAmount("Gypsum"))
        nAnhydrite = float(state.speciesAmount("Anhydrite"))
        mCa        = float(aprops.elementMolality("Ca"))

        print(f"{T:4.0f}  {nGypsum:8.4e} {nAnhydrite:10.4e} {mCa:8.4e}")

    return nGypsum, nAnhydrite, mCa

# Define the range of temperatures and pressure for the equilibrium calculations
T = np.arange(25.0, 76.0, 15.0)

# Output species amount after equilibration for a range of the
print(" -------------------------------------------")
print("  Final species amounts w.r.t. temperatures")
print(" -------------------------------------------")
print("   T   d(Gypsum) d(Anhydrite)     mCa")

# Run simmulations for P = 1 atm
P1 = 1.0
results = [minerals_in_purewater(x, P1) for x in T]  # [0] is needed to get the value of autodiff.real
nGypsumP1    = [molals[0] for molals in results]
nAnhydriteP1 = [molals[1] for molals in results]
mCaP1        = [molals[2] for molals in results]

# Run simmulations for P = 100 atm
P100 = 100.0
results      = [minerals_in_purewater(x, P100) for x in T]  # [0] is needed to get the value of autodiff.real
nGypsumP100    = [molals[0] for molals in results]
nAnhydriteP100 = [molals[1] for molals in results]
mCaP100        = [molals[2] for molals in results]

# Run simmulations for P = 1e3 atm
P1000 = 1000.0
results = [minerals_in_purewater(x, P1000) for x in T]  # [0] is needed to get the value of autodiff.real
nGypsumP1000    = [molals[0] for molals in results]
nAnhydriteP1000 = [molals[1] for molals in results]
mCaP1000       = [molals[2] for molals in results]

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
ax.plot(T, mCaP1000, label='P = 1000 amt', color="C3")
ax.plot(T, mCaP100, label='P = 100 amt', color="C2")
ax.plot(T, mCaP1, label='P = 1 amt', color="C1")
ax.legend(loc="best")
ax.set_title('Solubility of Gypsum and Anhydrite')
ax.grid()
ax.set_ylabel('Molality of Ca [molal]')
fig.savefig('solubility-gypsum-anhydrite-with-solver.png')
fig.tight_layout()

fig, ax = plt.subplots(1, 1)
ax.plot(T, nGypsumP1, label='Gypsum, P = 1 atm', color="C2", linestyle=":")
ax.plot(T, nAnhydriteP1, label='Anhydrite, P = 1 atm', color="C2", linestyle="-")
ax.plot(T, nGypsumP100, label='Gypsum, P = 1000 atm', color="C4", linestyle=":")
ax.plot(T, nAnhydriteP100, label='Anhydrite, P = 1000 atm', color="C4", linestyle="-")
ax.legend(loc="best")
ax.set_title('Gypsum and Anhydrite')
ax.grid()
ax.set_ylabel('Amounts of minerals [moles]')
fig.savefig('gypsum-anhydrite-with-solver.png')
fig.tight_layout()
