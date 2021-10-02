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
solution = AqueousPhase(speciate("H O C Ca Na Mg Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define an ion exchange phase
exchange_species = db.species().withAggregateState(AggregateState.IonExchange)
exchange = IonExchangePhase(speciesListToStringList(exchange_species))
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
state.setSpeciesAmount("Na+" , 1.10, "mol")
state.setSpeciesAmount("Mg+2", 0.48, "mol")
state.setSpeciesAmount("Ca+2", 1.90, "mol")
state.setSpeciesAmount("X-"  , 0.06, "mol")

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(specs)
solver.solve(state, conditions)
state.output("state-python-species-full.txt")

# Compute chemical properties for the aqueous phase at equilibrium state and output them
aprops = AqueousProps(state)
aprops.output("aprops-python-species-full.txt")
