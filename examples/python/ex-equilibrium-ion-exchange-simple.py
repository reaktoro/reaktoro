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
#   • Svetlana Kyas (29 September)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

def speciesListToStringList(specieslist):
    species_list = ""
    for s in specieslist:
        species_list += s.name() + " "
    return species_list

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2")
solution.setActivityModel(chain(
    ActivityModelHKF()
))

# Define an ion exchange phase
exchange_species = "NaX KX CaX2 MgX2"
exchange = IonExchangePhase(exchange_species)
exchange.setActivityModel(ActivityModelIonExchangeGainesThomas())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

# Specify conditions to be satisfied at the chemical equilibrium
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define conditions to be satisfied at chemical equilibrium
conditions = EquilibriumConditions(specs)
conditions.temperature(T, "celsius")
conditions.pressure(P, "bar")

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "bar")
state.setSpeciesMass("H2O"   , 1.00, "kg")
state.setSpeciesAmount("Na+" , 1.00, "mol")
state.setSpeciesAmount("K+"  , 1.00, "mol")
state.setSpeciesAmount("Mg+2", 1.00, "mol")
state.setSpeciesAmount("Ca+2", 1.00, "mol")
state.setSpeciesAmount("NaX" , 0.06, "umol") # set small to make sure we have plenty of water for available exchanger X-

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(specs)
solver.solve(state, conditions)
print(state)

aqprops = AqueousProps(state)
print("I  = %f mol/kgw" % aqprops.ionicStrength()[0])
print("pH = %f" % aqprops.pH()[0])
print("pE = %f" % aqprops.pE()[0])

exprops = IonExchangeProps(state)
print(exprops)