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
#   • Allan Leal (21 September 2022)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Note: ODML is currently implemented with only an on-demand clustering
# strategy for storing and searching learning data. However, this is not
# suitable for the example below, which would possibly benefit more from
# a nearest-neighbor search algorithm instead, to be implemented.
# -----------------------------------------------------------------------------

from reaktoro import *
from reaktplot import *
import numpy as np

# -----------------------------------------------------------------------------
# CONFIGURING CHEMICAL SYSTEM FOR THE CALCULATIONS
# -----------------------------------------------------------------------------

db = SupcrtDatabase("supcrtbl")

solution = AqueousPhase("H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-2 CO2(aq)")
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

gases = GaseousPhase("CO2(g) H2O(g)")
gases.setActivityModel(ActivityModelPengRobinson())

system = ChemicalSystem(db, solution, gases)

# -----------------------------------------------------------------------------
# PARAMETERS FOR THE CALCULATIONS
# -----------------------------------------------------------------------------

N = 50  # number of temperature and pressure points

temperatures = np.linspace(25.0, 90.0, N)  # temperatures from 25 to 90 °C
pressures = np.linspace(1.0, 300.0, N)  # pressures from 1 to 300 bar

# -----------------------------------------------------------------------------
# CONFIGURATION OF THE SMART AND EXACT EQUILIBRIUM SOLVERS
# -----------------------------------------------------------------------------

smart_solver = SmartEquilibriumSolver(system)
exact_solver = EquilibriumSolver(system)

exact_solver_options = EquilibriumOptions()

smart_solver_options = SmartEquilibriumOptions()
smart_solver_options.reltol = 0.001

exact_solver.setOptions(exact_solver_options)
smart_solver.setOptions(smart_solver_options)

# The initial chemical state in disequilibrium
state0 = ChemicalState(system)
state0.set("H2O(aq)", 1.0, "kg")
state0.set("Na+",     1.0, "mol")
state0.set("Cl-",     1.0, "mol")
state0.set("CO2(g)", 10.0, "mol")

exact_mCO2 = np.zeros((N, N))
smart_mCO2 = np.zeros((N, N))

xpredicted, ypredicted = [], []  # the points (x, y) = (T, P) where prediction occurred
xlearning, ylearning = [], []  # the points (x, y) = (T, P) where learning occurred

# -----------------------------------------------------------------------------
# CALCULATIONS USING SMART SOLVER RELYING ON ON-DEMAND MACHINE LEARNING (ODML)
# -----------------------------------------------------------------------------

print("Solving equilibrium problems using EquilibriumSolverODML...")

stopwatch = Stopwatch()

state = ChemicalState(state0)

for i, T in enumerate(temperatures):
    for j, P in enumerate(pressures):

        state.temperature(T, "celsius")
        state.pressure(P, "bar")

        stopwatch.start()

        result = smart_solver.solve(state)

        stopwatch.pause()

        assert result.succeeded(), f"Equilibrium calculation (smart) failed at {T} celsius and {P} bar."

        smart_mCO2[j, i] = state.props().elementAmountInPhase("C", "AqueousPhase")

        if result.predicted():
            xpredicted.append(T)
            ypredicted.append(P)
        else:
            xlearning.append(T)
            ylearning.append(P)

computing_time_smart_solver = stopwatch.time()

print(f"Finished! Computing time using SmartEquilibriumSolver: {computing_time_smart_solver} s")

# -----------------------------------------------------------------------------
# CALCULATIONS USING CONVENTIONAL/EXACT SOLVER
# -----------------------------------------------------------------------------

print("Solving equilibrium problems using EquilibriumSolver...")

stopwatch = Stopwatch()

state = ChemicalState(state0)

for i, T in enumerate(temperatures):
    for j, P in enumerate(pressures):

        state.temperature(T, "celsius")
        state.pressure(P, "bar")

        stopwatch.start()

        result = exact_solver.solve(state)

        stopwatch.pause()

        assert result.succeeded(), f"Equilibrium calculation (conventional) failed at {T} celsius and {P} bar."

        exact_mCO2[j, i] = state.props().elementAmountInPhase("C", "AqueousPhase")

computing_time_exact_solver = stopwatch.time()

print(f"Finished! Computing time using EquilibriumSolver: {computing_time_exact_solver} s")

# -----------------------------------------------------------------------------
# PLOTTING THE RESULTS
# -----------------------------------------------------------------------------

speedup = computing_time_exact_solver / computing_time_smart_solver
prediction_success_rate = len(xpredicted)/(N*N) * 100.0  # in percent

fig = Figure()
fig.title(f"CO2 SOLUBILITY IN BRINE - WITH ODML - SPEEDUP: {speedup:.1f} \n PREDICTED: {prediction_success_rate:.1f}%")
fig.xaxisTitle("Temperature [°C]")
fig.yaxisTitle("Pressure [bar]")
fig.drawContour(temperatures, pressures, smart_mCO2, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)))
fig.drawMarkers(xpredicted, ypredicted, "Prediction", MarkerSpecs().color("black").size(5))
fig.show()
fig.save("ex-smart-equilibrium-co2-solubility-brine-odml-yes.pdf")

fig = Figure()
fig.title("CO2 SOLUBILITY IN BRINE - WITHOUT ODML")
fig.xaxisTitle("Temperature [°C]")
fig.yaxisTitle("Pressure [bar]")
fig.drawContour(temperatures, pressures, exact_mCO2, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)))
fig.show()
fig.save("ex-smart-equilibrium-co2-solubility-brine-odml-no.pdf")

error = np.abs(exact_mCO2 - smart_mCO2)/np.abs(exact_mCO2) * 100.0  # error in percent

fig = Figure()
fig.title("CO2 SOLUBILITY IN BRINE - ERROR USING ODML")
fig.xaxisTitle("Temperature [°C]")
fig.yaxisTitle("Pressure [bar]")
fig.drawContour(temperatures, pressures, error, ContourSpecs().coloringModeHeatmap().line(LineSpecs().width(0)))
fig.show()
fig.save("ex-smart-equilibrium-co2-solubility-brine-odml-error.pdf")
