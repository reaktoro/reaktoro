# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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
#   • Svetlana Kyas (23 November, 2021)
#
# and since revised by:
#   • Allan Leal (28 August 2023)
#     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
# -----------------------------------------------------------------------------


from reaktoro import *
import numpy as np
from pathlib import Path

filepath = Path(__file__).parent.parent/"resources/phreeqc-rk-isotherm.dat"
db = PhreeqcDatabase.fromFile(str(filepath))

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2")
solution.set(ActivityModelPhreeqc(db))

# Define an ion exchange phase
exchange = IonExchangePhase("NaX KX CaX2")
exchange.set(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in atm

# Define equilibrium solver
solver = EquilibriumSolver(system)

# Define aqueous and exchange properties
exprops = IonExchangeProps(system)
aqprops = AqueousProps(system)

mols_K  = np.flip(np.linspace(0, 0.5, num=21))
mols_Ca = np.linspace(0, 0.25, num=21)

def equilibrate(solver, m_K, m_Ca):

    # Define initial equilibrium state
    state = ChemicalState(system)
    state.temperature(T, "celsius")
    state.pressure(P, "atm")
    state.set("H2O"   , 1.0 , "kg")
    state.set("K+"  , m_K , "mol")
    state.set("Ca+2", m_Ca, "mol")
    # Exchanger site
    state.set("NaX", 0.5, "mol")

    # Equilibrate chemical state
    solver.solve(state)

    # Update exchange and aqueous properties
    exprops.update(state)
    aqprops.update(state)

    return [float(state.speciesAmount("K+")),
            float(state.speciesAmount("Ca+2")),
            float(state.speciesAmount("KX")),
            float(state.speciesAmount("CaX2")),
            float(state.speciesAmount("NaX"))]

print(f"  input K+ input Ca+2      mol_K     mol_Ca     mol_KX   mol_CaX2    mol_NaX")

for mol_K, mol_Ca in zip(mols_K, mols_Ca):

    m_K, m_Ca, m_KX, m_CaX2, m_NaX = equilibrate(solver, mol_K, mol_Ca)
    print(f"{mol_K:6.4e} {mol_Ca:6.4e} {m_K:6.4e} {m_Ca:6.4e} {m_KX:6.4e} {m_CaX2:6.4e} {m_NaX:6.4e}")
