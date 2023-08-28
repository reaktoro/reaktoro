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
#   • Svetlana Kyas (23 November 2021)
#
# and since revised by:
#   • Allan Leal (28 August 2023)
#     - Using ActivityModelPhreeqc instead of ActivityModelHKF for aqueous phase.
# -----------------------------------------------------------------------------


import numpy as np
from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2")
solution.set(ActivityModelPhreeqc(db))

# Define an ion exchange phase
exchange = IonExchangePhase("NaX KX CaX2")
exchange.set(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define equilibrium solver
solver = EquilibriumSolver(system)

# Define aqueous and exchange properties
chemprops = ChemicalProps(system)
exprops = IonExchangeProps(system)

def equilibrate_state(m_Ca):

    # Define an equilibrium state
    state = ChemicalState(system)
    state.temperature(T, "celsius")
    state.pressure(P, "bar")
    state.set("H2O"   , 1.0 , "kg")
    state.set("K+"  , 0.1 , "mol")
    state.set("Ca+2", m_Ca, "mol")
    # Exchanger site
    state.set("NaX"  , 0.417, "mol")

    # Equilibrate chemical state
    solver.solve(state)

    # Update exchange and aqueous properties
    chemprops.update(state)
    exprops.update(state)

    print(f"{float(chemprops.elementAmount('K')):6.4e} "
          f"{float(chemprops.elementAmount('Ca')):6.4e} "
          f"{float(state.speciesAmount('K+')):6.4e} "
          f"{float(state.speciesAmount('Ca+2')):6.4e} "
          f"{float(state.speciesAmount('KX')):6.4e} "
          f"{float(state.speciesAmount('CaX2')):6.4e} "
          f"{float(state.speciesAmount('NaX')):6.4e}")

m_Ca = np.flip(np.linspace(0.0, 0.05, 21))

print("   input K   input Ca         K+       Ca+2         KX       CaX2        NaX")

for m in m_Ca:
    equilibrate_state(m)
