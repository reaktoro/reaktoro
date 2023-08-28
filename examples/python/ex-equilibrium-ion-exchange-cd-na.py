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
#   • Svetlana Kyas (23 July 2021)
#
# and since revised by:
#   • Allan Leal (28 August 2023)
#     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
# -----------------------------------------------------------------------------


from reaktoro import *
import numpy as np
import math

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O C Ca Na Cl Cd"))
solution.set(ActivityModelPhreeqc(db))

# Define an ion exchange phase
exchange = IonExchangePhase("NaX CdX2")
exchange.set(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in atm

# metal = {"Pb": "Pb+2",
#          "Cd": "Cd+2",
#          "Zn": "Zn+2",
#          "Cu": "Cu+2",
#          "Sr": "Sr+2",
#          "Fe": "Fe+3",
#          "Ba": "Ba+2",
#          "Ca": "Ca+2"}

metal = {"Cd": "Cd+2"}
metal_b = "Cd"

# Define equilibrium solver
solver = EquilibriumSolver(system)

# Define aqueous and exchange properties
props = ChemicalProps(system)

nums = 21
mols_NaX = np.linspace(0.01, 10.0, num=nums)
kds = np.zeros(nums)

def equilibrate(m_NaX):

    # Define initial equilibrium state
    state = ChemicalState(system)
    state.temperature(T, "celsius")
    state.pressure(P, "atm")
    state.set("H2O", 1.0 , "kg")
    state.set(metal[metal_b], 1.0, "mol")
    # Exchanger site
    state.set("NaX", m_NaX, "mol")

    # Equilibrate chemical state
    solver.solve(state)

    # Update exchange and aqueous properties
    props.update(state)

    # Fetch amount of sorbed and dissolved mineral
    b_tot  = float(props.elementAmount(metal_b))
    b_aq   = float(props.elementAmountInPhase(metal_b, "AqueousPhase"))
    b_surf = float(props.elementAmountInPhase(metal_b, "IonExchangePhase"))
    if b_surf != 0.0:
        kd = b_surf/b_aq
    else:
        kd = math.nan

    return [m_NaX, state.speciesAmount("NaX")[0], state.speciesAmount("CdX2")[0], kd]

print(f"input(NaX)     m(NaX)    m(CdX2) kd = m(CdX2)/m(Cd+2)")
for i, mol_NaX in enumerate(mols_NaX):

    m0_NaX, m_NaX, m_CdX2, kd = equilibrate(mol_NaX)
    print(f"{m0_NaX:6.4e} {m_NaX:6.4e} {m_CdX2:6.4e} {kd:20.4e}")

    # Save the distribution coefficient
    kds[i] = kd

from matplotlib import pyplot as plt
plt.plot(mols_NaX, kds)
plt.title("Distribution coefficient")
plt.xlabel('Initial amount of NaX [mol]')
plt.ylabel(r'K$_d$ = s$_I$ / c$_I$ [-]')
plt.grid()
plt.savefig("kd-cd-vs-initial-nax.png", bbox_inches='tight')
plt.close()
