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
#   • Allan Leal (16 July 2021)
#
# and since revised by:
#   • Allan Leal (22 July 2021)
# -----------------------------------------------------------------------------


from reaktoro import *


db = SupcrtDatabase("supcrtbl")

solution = AqueousPhase(speciate("Na Cl C Ca Mg Si"), exclude("organic"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

minerals = MineralPhases("Halite Calcite Magnesite Dolomite Quartz")

system = ChemicalSystem(db, solution, gases, minerals)

statex = ChemicalState(system)
statex.temperature(60.0, "celsius")
statex.pressure(100.0, "bar")
statex.set("H2O(aq)"  , 1.00, "kg")
statex.set("Halite"   , 1.00, "mol")
statex.set("Calcite"  , 1.00, "mol")
statex.set("Magnesite", 1.00, "mol")
statex.set("Quartz"   , 1.00, "mol")

equilibrate(statex)

propsx = ChemicalProps(statex)
Vx = propsx.volume()
Ux = propsx.internalEnergy()

statex.output("state-expected.txt")
propsx.output("props-expected.txt")

specs = EquilibriumSpecs(system)

idxV = specs.addInput("V")
idxU = specs.addInput("U")

volumeConstraint = ConstraintEquation()
volumeConstraint.id = "VolumeConstraint"
volumeConstraint.fn = lambda state, w: state.props().volume() - w[idxV]

internalEnergyConstraint = ConstraintEquation()
internalEnergyConstraint.id = "InternalEnergyConstraint"
internalEnergyConstraint.fn = lambda state, w: state.props().internalEnergy() - w[idxU]

specs.addConstraint(volumeConstraint)
specs.addConstraint(internalEnergyConstraint)

conditions = EquilibriumConditions(specs)
conditions.set("V", Vx)
conditions.set("U", Ux)
conditions.setLowerBoundPressure(1.0, "bar")

state = ChemicalState(system)
state.temperature(25.0, "celsius")
state.pressure(1.0, "bar")
state.set("H2O(aq)"  , 1.00, "kg")
state.set("Halite"   , 1.00, "mol")
state.set("Calcite"  , 1.00, "mol")
state.set("Magnesite", 1.00, "mol")
state.set("Quartz"   , 1.00, "mol")

solver = EquilibriumSolver(specs)

solver.solve(state, conditions)

props = ChemicalProps(state)

state.output("state.txt")
props.output("props.txt")

print("Success! Check outputted files `state.txt`, `props.txt`, `state-expected.txt`, `props-expected.txt`.")
print("Verify if `props.txt` and `props-expected.txt` are numerically equivalent.")
