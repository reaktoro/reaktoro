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
#   • Allan Leal (19 July 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = ThermoFunDatabase("aq17")

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2 H2O")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases()

system = ChemicalSystem(db, solution, gases, minerals)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

solver = EquilibriumSolver(specs)

conditions = EquilibriumConditions(specs)
conditions.temperature(60.0, "celsius")
conditions.pressure(100.0, "bar")
conditions.pH(4.0)
conditions.startWith("H2O@", 1.0, "kg")
conditions.startWith("Na+",  1.0, "mol")
conditions.startWith("Cl-",  1.0, "mol")
conditions.startWith("CO2", 10.0, "mol")

state = ChemicalState(system)

solver.solve(state, conditions)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")
