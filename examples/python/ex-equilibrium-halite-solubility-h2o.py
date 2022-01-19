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
#   • Svetlana Kyas (14 July 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

# Initialize a thermodynamic database
db = Database("supcrt98.yaml")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O Na Cl"))
aqueousphase.setActivityModel(ActivityModelHKF())

# Create a gaseous phase
mineralphase = MineralPhase("Halite")

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(mineralphase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Initial amount of NaCl
m0Halite = 100.0

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(25.0, "celsius")
state.setPressure(1.0, "bar")
state.setSpeciesMass("H2O(aq)", 1.0, "kg")
state.setSpeciesAmount("Halite", m0Halite, "mol")

# Define equilibrium solver and equilibrate given initial state
solver = EquilibriumSolver(system)
solver.solve(state)

# Output chemical state into the txt-file
state.output("state.txt")

print(f"Solubility of Halite in water is {(m0Halite - state.speciesAmount('Halite')) / state.speciesMass('H2O(aq)')} molal")
