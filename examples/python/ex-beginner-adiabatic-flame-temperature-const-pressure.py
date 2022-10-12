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
#   • Allan Leal (19 July 2021)
# -----------------------------------------------------------------------------

from reaktoro import *


reactantName = "CH4"
oxidizerName = "Air"

reactantAmount = 1.0 # in mol
oxidizerAmount = 2.0 # in mol

print("LOADING NASA THERMO.INP DATABASE...")
db = NasaDatabase("nasa-cea")

reactant = db.species().getWithName(reactantName)
oxidizer = db.species().getWithName(oxidizerName)

symbols = list(set(reactant.elements().symbols() + oxidizer.elements().symbols()))

gases = GaseousPhase(speciate(symbols))
condensed = CondensedPhases(speciate(symbols))

print("CREATING CHEMICAL SYSTEM...")

system = ChemicalSystem(db, gases, condensed)

print(f"CHEMICAL SYSTEM CREATED WITH {system.elements().size()} ELEMENTS AND {system.species().size()} SPECIES")

state0 = ChemicalState(system)
state0.temperature(25.0, "celsius")
state0.pressure(1.0, "atm")
state0.setSpeciesAmounts(1e-16)
state0.set(reactantName, reactantAmount, "mol")
state0.set(oxidizerName, oxidizerAmount, "mol")

props0 = ChemicalProps(state0)

specs = EquilibriumSpecs(system)
specs.pressure()
specs.enthalpy()

conditions = EquilibriumConditions(specs)
conditions.pressure(props0.pressure())
conditions.enthalpy(props0.enthalpy())
conditions.setLowerBoundTemperature(298.15, "celsius")
conditions.setUpperBoundTemperature(4000.0, "celsius")

state = ChemicalState(state0)

print("PERFORMING THE CHEMICAL EQUILIBRIUM CALCULATION...")
solver = EquilibriumSolver(specs)

result = solver.solve(state, conditions)

assert result.succeeded(), "The calculation did not succeed!"

print(state)
print(f"COMPUTED ADIABATIC FLAME TEMPERATURE: {state.temperature()} K")
print(f"ITERATIONS: {result.iterations()}")
