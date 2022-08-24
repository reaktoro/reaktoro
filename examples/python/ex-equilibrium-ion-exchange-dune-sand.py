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
#   • Svetlana Kyas (14 October 2021)
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O C Ca Na Mg Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define an ion exchange phase
exchange_species = "NaX CaX2 MgX2"
exchange = IonExchangePhase(exchange_species)
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define initial equilibrium state
state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
# Scale solution recipe to match the values of the PHREEQC examples
state.set("H2O"   , 1.e6, "kg")
state.set("Na+" , 1.10, "kmol")
state.set("Mg+2", 0.48, "kmol")
state.set("Ca+2", 1.90, "kmol")
# Set the number of exchange assuming that it is completely occupied by sodium
state.set("NaX" , 0.06, "mol")

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)
solver.solve(state)
print(state)

aqprops = AqueousProps(state)
print("I  = %f mol/kgw" % aqprops.ionicStrength()[0])
print("pH = %f" % aqprops.pH()[0])
print("pE = %f" % aqprops.pE()[0])

chemprops = ChemicalProps(state)
exprops = IonExchangeProps(chemprops)
print(exprops)
