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
import numpy as np

# Function to calculate solubility of CO2 in the NaCl-brine
def solubility_co2(system, T, P, mNaCl):

    # Initial amount of CO2 gas
    n0CO2g = 10.0

    # Define initial chemical state corresponding to the NaCl-brine of the given concentration
    state = ChemicalState(system)
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesMass("H2O", 1.0, "kg")
    state.setSpeciesAmount("CO2(g)", n0CO2g, "mol")
    state.setSpeciesAmount("Na+", mNaCl, "mol")
    state.setSpeciesAmount("Cl-", mNaCl, "mol")

    # Calculate equilibrium state
    solver = EquilibriumSolver(system)
    res = solver.solve(state)

    # Throw exception if the equilibrium couldn't be found
    if not res.optima.succeeded:
        raise RuntimeError("Equilibrium calculation did not succeed!")

    # Fetch resulting aqueous properties of the chemical state
    aqprops = AqueousProps(state)

    # Return concentration of the carbon in the aqueous phase
    return aqprops.elementMolality("C")[0]

# Initialize a thermodynamic database
db = PhreeqcDatabase("phreeqc.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
aqueousphase.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2"),
))

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

# Collecting all above-defined phases
phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)

# Construct the chemical system
system = ChemicalSystem(phases)

# Define the range of temperatures and pressure for the equilibrium calculations
T = np.arange(25.0, 90.0, 5.0)
P = 100.0

# Calculate CO2 solubilities for the range of the temperatures and different brine concentrations
mCO2_1 = [solubility_co2(system, x, P, mNaCl=1.0) for x in T]
mCO2_2 = [solubility_co2(system, x, P, mNaCl=2.0) for x in T]
mCO2_3 = [solubility_co2(system, x, P, mNaCl=4.0) for x in T]

# Output the results
print("# ---------------------------------------------------------------")
print("# CO2 solubilities w.r.t. temperatures")
print("# ---------------------------------------------------------------")
print(" T      1 mol NaCl-brine   2 mol NaCl-brine   4 mol NaCl-brine")
for i in range(len(T)):
    print(f"{T[i]:3.0f} : {mCO2_1[i]:18.4f} {mCO2_2[i]:18.4f} {mCO2_3[i]:18.4f}")
