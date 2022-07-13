# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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
#   • Svetlana Kyas (13 July 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *
import numpy as np

# Load database
db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O C Na N Cl Ca Mg"))
solution.setActivityModel(ActivityModelHKF())

# Define an ion exchange phase
exchange = IonExchangePhase("NaX CaX2")
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Define an ion exchange phase
mineral = MineralPhase("Calcite")

# Create chemical system
system = ChemicalSystem(db, solution, exchange, mineral)

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)

# Ion-exchange and aqueous properties
exprops = IonExchangeProps(system)
aqprops = AqueousProps(system)

# Auxiliary arrays and constants
num_ph = 31
num_exchanges = 4
mols_NH4 = np.linspace(0, 60.0, num=num_ph)
mols_NaX = 1e-3 * np.array([0, 2.5, 12.5, 25])

pHs          = np.zeros((num_exchanges, num_ph))
Is           = np.zeros((num_exchanges, num_ph))
mols_Ca      = np.zeros((num_exchanges, num_ph))
mols_CaX2    = np.zeros((num_exchanges, num_ph))
mols_Calcite = np.zeros((num_exchanges, num_ph))

def equilibrate(solver, mol_NaX, mol_NH4):

    # Initial values of the carbonates and mass of water
    m0Calcite = 10.0

    # Define initial equilibrium state
    state = ChemicalState(system)
    state.setTemperature(25.0, "celsius")
    state.setPressure(1, "atm")
    # Seawater
    state.set("H2O" , 1.0, "kg")
    state.set("Na+" , 1.10, "mmol")
    state.set("Mg+2", 0.48, "mmol")
    state.set("Ca+2", 1.90, "mmol")
    state.set("NH4+", mol_NH4, "mmol")
    # Ammonia
    state.set("Calcite", m0Calcite, "mol")

    # Calculate chemical state corresponding to the seawater
    solver.solve(state)

    # Update aqueous properties
    aqprops.update(state)
    pH = aqprops.pH()[0]

    # Exchanger
    state.setSpeciesAmount("NaX" , mol_NaX, "mol")

    # Equilibrate the seawater with carbonates
    solver.solve(state)

    # Update xchange and aqueous properties
    exprops.update(state)
    aqprops.update(state)

    mCalcite = state.speciesAmount("Calcite")

    return [mol_NH4,
            pH,
            aqprops.ionicStrength()[0],
            state.speciesAmount("Ca+2")[0],
            m0Calcite - mCalcite[0],
            state.speciesAmount("CaX2")[0]]

exchange_counter = 0
for mol_NaX in mols_NaX:

    print(f"----------------------------------------------------------\n"
          f"mol_NaX = {mol_NaX:6.4e}\n"
          f"----------------------------------------------------------\n"
          f" i  mol_NH4   pH        I       m_Ca  m_Calcite     m_CaX2  ")
    ph_counter = 0
    for mol_NH4 in mols_NH4:

        mol_NH4, pH, I, m_Ca, m_Calcite, m_CaX2 = equilibrate(solver, mol_NaX, mol_NH4)
        pHs[exchange_counter, ph_counter]          = pH
        Is[exchange_counter, ph_counter]           = I
        mols_Ca[exchange_counter, ph_counter]      = m_Ca
        mols_Calcite[exchange_counter, ph_counter] = m_Calcite
        mols_CaX2[exchange_counter, ph_counter]    = m_CaX2

        print(f"{ph_counter:2d} {mol_NH4:8.2f} {pH:4.2f} {I:6.2e} {m_Ca:6.4e} {m_Calcite:6.4e} {m_CaX2:6.4e} ")

        ph_counter += 1

    exchange_counter += 1


import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
for i in range(0, num_exchanges):
    plt.plot(pHs[i, :], mols_Calcite[i, :], label=f'{mols_NaX[i]} mol of NaX', color=colors[i])
plt.legend(loc="best")
plt.xlabel('pH [-]')
plt.ylabel('Amount [mol]')
plt.title('Amount of dissolved Calcite')
plt.grid()
plt.savefig('dissolved-calcite-vs-ph-NaX.png', bbox_inches='tight')
plt.close('all')

plt.figure()
for i in range(0, num_exchanges):
    plt.plot(pHs[i, :], mols_Ca[i, :], label=f'{mols_NaX[i]} mol of NaX', color=colors[i])
plt.legend(loc="best")
plt.xlabel('pH [-]')
plt.ylabel('Amount [mol]')
plt.title('Amount of Ca+2')
plt.grid()
plt.savefig('dissolved-ca-vs-ph-NaX.png', bbox_inches='tight')
plt.close()

plt.figure()
for i in range(0, num_exchanges):
    plt.plot(pHs[i, :], mols_CaX2[i, :], label=f'{mols_NaX[i]} mol of NaX', color=colors[i])
plt.legend(loc="best")
plt.xlabel('pH [-]')
plt.ylabel('Amount [mol]')
plt.title('Amount of CaX2')
plt.grid()
plt.savefig('dissolved-cax2-vs-ph-NaX.png', bbox_inches='tight')
plt.close('all')
