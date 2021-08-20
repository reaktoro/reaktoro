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
#   • Allan Leal (4 August 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------

from reaktoro import *


db = SupcrtDatabase("supcrtbl")

solution = AqueousPhase(speciate("Na Cl C"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

system = ChemicalSystem(db, solution, gases)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.phaseAmount("GaseousPhase")

solver = EquilibriumSolver(specs)

state = ChemicalState(system)
state.setTemperature(50.0, "celsius")
state.setPressure(300.0, "bar")
state.set("H2O(aq)", 1.0, "kg")
state.set("Na+",     1.0, "mol")
state.set("Cl-",     1.0, "mol")
state.set("CO2(aq)", 1.0, "mol")

conditions = EquilibriumConditions(specs)
conditions.temperature(50.0, "celsius")
conditions.phaseAmount("GaseousPhase", 1e-10, "mol")
conditions.setLowerBoundPressure(1.0, "bar")
conditions.setUpperBoundPressure(1000.0, "bar")

solver.solve(state, conditions)

print(f"Pressure: {state.pressure() * 1e-5} bar")
