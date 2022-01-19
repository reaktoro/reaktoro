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
#   • Svetlana Kyas (29 September)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g)")
gases.setActivityModel(ActivityModelPengRobinson())

system = ChemicalSystem(db, solution, gases)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
state.setSpeciesMass("H2O"     , 1.0 , "kg")
state.setSpeciesAmount("CO2(g)", 10.0, "mol")
state.setSpeciesAmount("Na+"   , 4.00, "mol")
state.setSpeciesAmount("Cl-"   , 4.00, "mol")

solver = EquilibriumSolver(system)
solver.solve(state)

props = ChemicalProps(state)
props.output("props.txt")

aprops = AqueousProps(state)
aprops.output("aprops.txt")

print("Success! Check outputted files `props.txt` and `aprops.txt`.")
