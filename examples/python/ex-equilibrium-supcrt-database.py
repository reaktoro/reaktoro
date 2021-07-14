# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
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

# Import reaktoro python package
from reaktoro import *

# Initialize a thermodynamic database
db = SupcrtDatabase("supcrt98.xml")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
aqueousphase.setActivityModel(ActivityModelHKF())

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(25.0, "celsius")
state.setPressure(1.0, "bar")
state.setSpeciesMass("H2O(l)", 1.0, "kg")
state.setSpeciesAmount("CO2(g)", 10.0, "mol")
state.setSpeciesAmount("Na+", 1.0, "mol")
state.setSpeciesAmount("Cl-", 1.0, "mol")

# Define equilibrium solver and equilibrate given initial state
solver = EquilibriumSolver(system)
solver.solve(state)

# Obtain species composition from the equilibrated state
n = state.speciesAmounts()

for i in range(len(n)):
    print(f"{system.species(i).name():12s} {n[i][0]:6.4e}")
