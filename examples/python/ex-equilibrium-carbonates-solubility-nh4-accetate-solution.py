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
#   • Svetlana Kyas (13 July, 2022)
# -----------------------------------------------------------------------------


from reaktoro import *

# Function to calculate equilibrium of the carbonates and seawater
def carbonates_in_seawater(T):

    # Initial values of the carbonates and mass of water with nh4-accetate
    n0Calcite = 10.0
    water_kg = 1.0

    # Define initial equilibrium state corresponding to the seawater
    state_sw = ChemicalState(system)
    state_sw.temperature(T, "celsius")
    state_sw.pressure(1.0, "atm")
    state_sw.set("H2O", water_kg, "kg")
    state_sw.set("NH4+", 1, "mol")

    # Calculate chemical state corresponding to the seawater
    res = solver.solve(state_sw)

    # Update aqueous properties
    aqprops.update(state_sw)
    print("pH (before adding Calcite) = ", aqprops.pH())

    # Throw exception if the equilibrium couldn't be found
    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    # Add carbonates
    state_sw.set("Calcite", n0Calcite, "mol")

    # Equilibrate the seawater with carbonates
    solver.solve(state_sw)

    # Updated aqueous properties
    aqprops.update(state_sw)
    print("pH (after adding Calcite)  = ", aqprops.pH())

    # Fetch values of the specified species
    nCalcite = state_sw.speciesAmount("Calcite")
    nCa2 = state_sw.speciesAmount("Ca+2")
    pH = aqprops.pH()

    return (pH[0], nCalcite[0], nCa2[0])

# Initialize a thermodynamic database
db = PhreeqcDatabase("phreeqc.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Cl Na K Mg S Si N"))
aqueousphase.setActivityModel(ActivityModelPitzerHMW())

# Create carbonates phases
calcitephase = MineralPhase("Calcite")
dolomitephase = MineralPhase("Dolomite")

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(calcitephase)
phases.add(dolomitephase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Create the equilibrium solver
solver = EquilibriumSolver(system)

# Create the equilibrium solver
aqprops = AqueousProps(system)

# Define the range of temperatures and pressure for the equilibrium calculations
temperatures = [25.0] #np.arange(25.0, 91.0, 5.0)
# import numpy as np
# temperatures = np.arange(25.0, 91.0, 5.0)

# Fetch specific species amounts
species_amounts = [carbonates_in_seawater(x) for x in temperatures]  # [0] is needed to get the value of autodiff.real
pH       = [molals[0] for molals in species_amounts]
mCalcite = [molals[1] for molals in species_amounts]
mCa2     = [molals[2] for molals in species_amounts]

# Output species amount after equilibration for a range of the
print(" --------------------------------------------------------------------------------")
print("  Final species amounts w.r.t. temperatures")
print(" --------------------------------------------------------------------------------")
print("   T      pH      Calcite        Ca++")
for i in range(len(temperatures)):
    print(f"{temperatures[i]:4.0f}  {pH[i]:6.4f} {mCalcite[i]:12.4f} {mCa2[i]:12.4e}")
