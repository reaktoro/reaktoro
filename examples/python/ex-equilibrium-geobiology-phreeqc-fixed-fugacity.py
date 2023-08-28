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
#   • Svetlana Kyas (21 July, 2022)
# -----------------------------------------------------------------------------

# NOTE (by Allan Leal, 28 August 2023): This example requires the Pitzer
# parameters in the PHREEQC database below to be used, but this is not being
# done here yet. The activity model should be changed to ActivityModelPitzer.

from reaktoro import *
import numpy as np
import pandas as pd
from pathlib import Path

filepath = Path(__file__).parent.parent/"resources/phreeqc-toner-catling.dat"
db = PhreeqcDatabase.fromFile(str(filepath))

# Define aqueous phase
solution = AqueousPhase(speciate("H O C Na Cl P"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define mineral phases
minerals = MineralPhases("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O Halite Na2(HPO4):7H2O")

# Define chemical system
system = ChemicalSystem(db, solution, minerals)

# Define aqueous and chemical properties
props = ChemicalProps(system)
aprops = AqueousProps(system)

# Define equilibrium specifications
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

# Define conditions to be satisfied at the chemical equilibrium state
conditions = EquilibriumConditions(specs)

# Define equilibrium conditions
opts = EquilibriumOptions()
opts.epsilon = 1e-13

# Auxiliary arrays
num_log10pCO2s = 71
co2pressures = np.flip(np.linspace(-5.0, 2.0, num=num_log10pCO2s))
temperatures = np.array([0, 25, 50])

# Output dataframe
data = pd.DataFrame(columns=["T", "ppCO2", "pH", "mCO3", "mHCO3", "x", "amount_P"])

for T in temperatures:
    for log10pCO2 in co2pressures:

        conditions.temperature(T, "celsius")
        conditions.pressure(1.0, "atm")
        conditions.fugacity("CO2", 10 ** log10pCO2, "bar")

        state = ChemicalState(system)
        state.set("H2O"           ,   1.0, "kg")
        state.set("Nahcolite"     ,  10.0, "mol")
        state.set("Halite"        ,  10.0, "mol")
        state.set("Na2(HPO4):7H2O",  10.0, "mol")
        state.set("CO2"           , 100.0, "mol")

        solver = EquilibriumSolver(specs)
        solver.setOptions(opts)

        # Equilibrate the solution with the given initial chemical state and desired conditions at the equilibrium
        res = solver.solve(state, conditions)
        # Stop if the equilibration did not converge or failed
        if res.failed(): continue

        props.update(state)
        aprops.update(state)

        mCO3 = float(state.speciesAmount("CO3-2"))
        mHCO3 = float(state.speciesAmount("HCO3-"))
        x = 100 * 2 * mCO3 / (mHCO3 + 2 * mCO3)

        data.loc[len(data)] = [T, log10pCO2, float(aprops.pH()),
                               mCO3, mHCO3, x,
                               float(props.elementAmountInPhase("P", "AqueousPhase"))]

import matplotlib.pyplot as plt
colors = ['teal', 'darkred', 'indigo', 'coral', 'rosybrown', 'steelblue', 'seagreen', 'palevioletred', 'darkred', 'darkkhaki', 'cadetblue', 'indianred']

plt.figure()

data_T = data[data["T"] == 0] # fetch the columns with Pb+2
ax = data_T.plot(x="ppCO2", y="pH", color=colors[0], label=f"T = {0} °C")
for T, color in zip(temperatures[1:], colors[1:]):
    data_T = data[data["T"] == T]
    data_T.plot(x="ppCO2", y="pH", ax=ax, color=color, label=f"T = {T} °C")
ax.set_xlabel(r'$\log_{10}(\rm{pCO}_2)$ [-]')
ax.set_ylabel('pH [-]')
ax.legend(loc="best")
plt.grid()
plt.savefig(f'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()

data_T = data[data["T"] == 0] # fetch the columns with Pb+2
ax = data_T.plot(x="ppCO2", y="amount_P", color=colors[0], label="T = 0 °C")
for T, color in zip(temperatures[1:], colors[1:]):
    data_T = data[data["T"] == T]
    data_T.plot(x="ppCO2", y="amount_P", ax=ax, color=color, label=f"T = {T} °C")
ax.set_xlabel(r'$\log_{10}(\rm{pCO}_2)$ [-]')
ax.set_ylabel('P amount [mol]')
ax.legend(loc='best')
plt.grid()
plt.savefig(f'P-vs-ppCO2.png', bbox_inches='tight')
plt.close()

data_T = data[data["T"] == 25]
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'$\frac{2[\sf{CO}_3^{-2}]}{[\sf{HCO}_3^-] + 2[\sf{CO}_3^{-2}]}$ [%]')
ax1.set_ylabel('pH [-]', color=colors[3])
data_T.plot(x="x", y="pH", ax=ax1, color=colors[3], label='pH vs x')
ax1.legend(loc='upper center')
ax1.tick_params(axis='y', labelcolor=colors[3])

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set_ylabel(r'$\log_{10}(\sf{pCO}_2)$ [-]', color=colors[5])  # we already handled the x-label with ax1
data_T.plot(x="x", y="ppCO2", ax=ax2, color=colors[5], label=r'$\log_{10}$(pCO$_2$) vs x')
ax2.tick_params(axis='y', labelcolor=colors[5])
ax2.legend(loc="lower center")
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.grid()
plt.savefig(f'pH-ppCO2-vs-x.png', bbox_inches='tight')
plt.close()
