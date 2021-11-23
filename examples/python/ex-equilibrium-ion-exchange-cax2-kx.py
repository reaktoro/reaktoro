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
#   • Svetlana Kyas (23 November, 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------
import numpy as np
from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2")
solution.setActivityModel(ActivityModelHKF())

# Define an ion exchange phase
exchange = IonExchangePhase("NaX KX CaX2")
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

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
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesMass("H2O"   , 1.0 , "kg")
    state.setSpeciesAmount("K+"  , 0.1 , "mol")
    state.setSpeciesAmount("Ca+2", m_Ca, "mol")
    # Exchanger site
    state.setSpeciesAmount("NaX"  , 0.417, "mol")

    # Equilibrate chemical state
    solver.solve(state)

    # Update exchange and aqueous properties
    chemprops.update(state)
    exprops.update(state)

    print(f"{chemprops.elementAmount('K')[0]:6.4e} {chemprops.elementAmount('Ca')[0]:6.4e} "
          f"{state.speciesAmount('K+')[0]:6.4e} {state.speciesAmount('Ca+2')[0]:6.4e} "
          f"{state.speciesAmount('KX')[0]:6.4e} {state.speciesAmount('CaX2')[0]:6.4e} "
          f"{state.speciesAmount('NaX')[0]:6.4e}")

m_Ca = np.flip(np.linspace(0.0, 0.05, 21))

print("   input K   input Ca         K+       Ca+2         KX       CaX2        NaX")

for m in m_Ca:
    equilibrate_state(m)
