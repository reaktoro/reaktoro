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


from reaktoro import *
import numpy as np
import pandas as pd

# Load Phreeqc database
db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O C Na N Cl Ca Mg"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define an ion exchange phase
exchange = IonExchangePhase("NaX CaX2")
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Define an ion exchange phase
mineral = MineralPhase("Calcite")

# Create the chemical system
system = ChemicalSystem(db, solution, exchange, mineral)

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)

# Define ion exchange properties and aqueous properties
aqprops = AqueousProps(system)

# Sampling arrays of NH4 ions' amounts and NaX exchange species
num_ph = 31
num_exchangers = 4
mols_NH4 = np.linspace(0, 60.0, num=num_ph)
mols_NaX = 1e-3 * np.array([0.0, 2.5, 12.5, 25])

# Output dataframe
data = pd.DataFrame(columns=["amount_NaX", "amount_NH4", "pH", "I",
                             "amount_Ca", "amount_CaX2", "delta_Calcite"])

for mol_NaX in mols_NaX:
    for mol_NH4 in mols_NH4:

        # Initial amount of Calcite
        m0Calcite = 10.0

        # Define initial equilibrium state
        state = ChemicalState(system)
        state.setTemperature(25.0, "celsius")
        state.setPressure(1.0, "atm")
        # Seawater
        state.set("H2O" , 1.0    , "kg")
        state.set("Na+" , 1.10   , "mmol")
        state.set("Mg+2", 0.48   , "mmol")
        state.set("Ca+2", 1.90   , "mmol")
        state.set("NH4+", mol_NH4, "mmol")
        # Ammonia
        state.set("Calcite", m0Calcite, "mol")

        # Calculate chemical state corresponding to the seawater
        res = solver.solve(state)
        # Stop if the equilibration hasn't converged
        if not res.optima.succeeded: continue

        # Update aqueous properties and evaluate pH
        aqprops.update(state)
        pH = float(aqprops.pH())

        # Exchanger
        state.set("NaX", mol_NaX, "mol")

        # Equilibrate the seawater with carbonates
        res = solver.solve(state)
        if not res.optima.succeeded: continue

        # Update aqueous properties to evaluate ionic strength
        aqprops.update(state)
        chemprops = state.props()

        # Collect the value to be added to the dataframe in the following order
        # "amount_NaX", "amount_NH4", "amount_pH", "I", "amount_Ca", "amount_CaX2", "delta_Calcite"
        data.loc[len(data)] = [mol_NaX, mol_NH4, pH, float(aqprops.ionicStrength()),
                              float(chemprops.elementAmountInPhase("Ca", "AqueousPhase")), float(state.speciesAmount("CaX2")),
                              m0Calcite - float(state.speciesAmount("Calcite"))]

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
data_NaX = data[data['amount_NaX'] == 0]
ax = data_NaX.plot(x="pH", y="delta_Calcite", color=colors[0], label=f"m(NaX) = {mols_NaX[0]} mol")
for mol_NaX, color in zip(mols_NaX[1:], colors[1:]):
    data_NaX = data[data['amount_NaX'] == mol_NaX]
    data_NaX.plot(x="pH", y="delta_Calcite", ax=ax, color=color, label=f"m(NaX) = {mol_NaX} mol")
ax.legend(loc="best")
ax.set_xlabel('pH [-]')
ax.set_ylabel('Amount [mol]')
ax.set_title('Amount of dissolved Calcite')
plt.grid()
plt.savefig('dissolved-calcite-vs-ph-NaX.png', bbox_inches='tight')
plt.close('all')

plt.figure()
data_NaX = data[data['amount_NaX'] == 0]
ax = data_NaX.plot(x="pH", y="amount_Ca", color=colors[0], label=f"m(NaX) = {mols_NaX[0]} mol")
for mol_NaX, color in zip(mols_NaX[1:], colors[1:]):
    data_NaX = data[data['amount_NaX'] == mol_NaX]
    data_NaX.plot(x="pH", y="amount_Ca", ax=ax, color=color, label=f"m(NaX) = {mol_NaX} mol")
plt.legend(loc="best")
ax.set_xlabel('pH [-]')
ax.set_ylabel('Amount [mol]')
ax.set_title(r'Amount of Ca$^{+2}$')
plt.grid()
plt.savefig('dissolved-ca-vs-ph-NaX.png', bbox_inches='tight')
plt.close()

plt.figure()
data_NaX = data[data['amount_NaX'] == 0]
ax = data_NaX.plot(x="pH", y="amount_CaX2", color=colors[0], label=f"m(NaX) = {mols_NaX[0]} mol")
for mol_NaX, color in zip(mols_NaX[1:], colors[1:]):
    data_NaX = data[data['amount_NaX'] == mol_NaX]
    data_NaX.plot(x="pH", y="amount_CaX2", ax=ax, color=color, label=f"m(NaX) = {mol_NaX} mol")
plt.legend(loc="best")
plt.xlabel('pH [-]')
plt.ylabel('Amount [mol]')
plt.title(r'Amount of CaX$_2$')
plt.grid()
plt.savefig('cax2-vs-ph-NaX.png', bbox_inches='tight')
plt.close('all')
