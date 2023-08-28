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
#   • Svetlana Kyas (21 July, 2022)
# -----------------------------------------------------------------------------

# NOTE (by Allan Leal, 28 August 2023): This example requires the Pitzer
# parameters in the PHREEQC database below to be used, but this is not being
# done here yet. The activity model should be changed to ActivityModelPitzer.

from reaktoro import *
import numpy as np
from pathlib import Path

# Fluorapatite
#     Ca5(F)(PO4)3 = 5Ca+2 + F- + 3PO4-3
#     log_k     -59.6
#     -analytical_expression -1917.945184 0 87834.57783 631.9611081 0 0
# Hydroxylapatite
#     Ca5(OH)(PO4)3 = 5Ca+2 + OH- + 3PO4-3
#     log_k     -58.517
#     -analytical_expression -1.6657 -0.098215 -8219.41 0 0 0
filepath = Path(__file__).parent.parent/"resources/phreeqc-extended.dat"
db = PhreeqcDatabase.fromFile(str(filepath))

# Define the aqueous phase
solution = AqueousPhase(speciate(StringList("H O C Na Cl Ca P")))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

# Define minerals' phases
minerals = MineralPhases("Calcite Fluorapatite Hydroxylapatite")

# Define the chemical system
system = ChemicalSystem(db, solution, minerals)

# Define aqueous and chemical properties
aprops = AqueousProps(system)
props = ChemicalProps(system)

# Define equilibrium specifications
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

# Define equilibrium conditions
conditions = EquilibriumConditions(specs)
conditions.pressure(1.0, "atm")

solver = EquilibriumSolver(specs)
opts = EquilibriumOptions()
opts.epsilon = 1e-13
solver.setOptions(opts)

num_temperatures = 11
num_log10pCO2s = 11
temperatures =  np.linspace(0.0, 50.0, num=num_temperatures)
co2ppressures = np.linspace(-4.0, 0.0, num=num_log10pCO2s)

data_size = 2
data_pH = np.zeros((num_temperatures, num_log10pCO2s))
data_P = np.zeros((num_temperatures, num_log10pCO2s))

for i in range(0, num_temperatures):
    for j in range(0, num_log10pCO2s):

        print(f"{i}, {j}: {temperatures[i]}, {co2ppressures[j]}")

        conditions.temperature(temperatures[i], "celsius")
        conditions.pressure(1.0, "atm")
        conditions.fugacity("CO2", 10 ** co2ppressures[j], 'atm')

        state = ChemicalState(system)
        state.set("H2O"            ,   1.0, "kg")
        state.set("Calcite"        ,  10.0, "mol")
        state.set("Fluorapatite"   ,  10.0, "mol")
        state.set("Hydroxylapatite",  10.0, "mol")
        state.set("CO2"            , 100.0, "mol")

        solver.solve(state, conditions)

        aprops.update(state)
        props.update(state)
        data_pH[i, j] = float(aprops.pH())
        data_P[i, j] = float(props.elementAmountInPhase("P", "AqueousPhase"))

import matplotlib as ml
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

levels = MaxNLocator(nbins=15).tick_values(data_pH.min(), data_pH.max())
cmap = ml.colormaps['rainbow']
reversed_map = cmap.reversed() # reversing the original colormap using reversed() function
norm = BoundaryNorm(levels, ncolors=reversed_map.N, clip=True)

fig, ax = plt.subplots(1, 1)
plt.rcParams['axes.grid'] = False
im = plt.pcolormesh(co2ppressures, temperatures, data_pH, cmap=reversed_map, norm=norm, shading='auto')
fig.colorbar(im, ax=ax)
ax.set_title('pH [-]')
ax.set_ylabel(r'T [$^\circ$C]')
ax.set_xlabel(r'$\log_{10}(\sf{pCO}_2)$ [-]')
plt.savefig('pH-colormesh.png', bbox_inches='tight')
plt.close('all')

levels = MaxNLocator(nbins=15).tick_values(np.log10(data_P).min(), np.log10(data_P).max())
cmap = ml.colormaps['cividis']
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, ax = plt.subplots(1, 1)
im = plt.pcolormesh(co2ppressures, temperatures, np.log10(data_P), cmap=cmap, norm=norm,  shading='auto')
fig.colorbar(im, ax=ax)
ax.set_title(r'Amount $\log_{10}(P)$ [mol]')
ax.set_ylabel(r'T [$^\circ$C]')
ax.set_xlabel(r'$\log_{10}(\sf{pCO}_2)$ [-]')
plt.savefig('log10P-colormesh.png', bbox_inches='tight')
plt.close('all')
