# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2020 Allan Leal
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


from reaktoro import *

db = SupcrtDatabase("supcrt98.xml")

aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
aqueousphase.setActivityModel(ActivityModelHKF())

gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)

system = ChemicalSystem(phases)

state = ChemicalState(system)

state.setTemperature(25.0, "celsius")
state.setPressure(1.0, "bar")
state.setSpeciesMass("H2O(l)", 1.0, "kg")
state.setSpeciesAmount("CO2(g)", 10.0, "mol")
state.setSpeciesAmount("Na+", 1.0, "mol")
state.setSpeciesAmount("Cl-", 1.0, "mol")

options = EquilibriumOptions()
options.optima.output.active = True
# options.optima.kkt.method = Optima.SaddlePointMethod.Fullspace
# options.optima.kkt.method = Optima.SaddlePointMethod.Nullspace
# options.optima.kkt.method = Optima.SaddlePointMethod.Rangespace
# options.optima.linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor = 0.0  # zero to disable line search
# options.optima.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = 0.0

solver = EquilibriumSolver(system)
solver.setOptions(options)

solver.solve(state)

n = state.speciesAmounts()

for i in range(len(n)):
    print(f"{system.species(i).name()} {n[i]}")
