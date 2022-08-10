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
#   • Svetlana Kyas (10 August 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O Cl Ca Fe"))
solution.setActivityModel(ActivityModelHKF())

# Define surface complexation  species list
species_list = db.species().withAggregateState(AggregateState.Adsorbed)
species_str = ' '.join(extractNames(species_list))

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined strong site of the complexation surface
site_Hfo_s = ComplexationSurfaceSite()
site_Hfo_s.setName("Hfo_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite(site_Hfo_s)

# Defined weak site of the complexation surface
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(species_list)

# Add specified surface as parameters for the activity model for the complexation surface
params = ActivityModelSurfaceComplexationParams()
params.surface = surface_Hfo

# Define surface complexation phase and set an activity model
complexation_phase_Hfo = SurfaceComplexationPhase(species_str)
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params))

# Create chemical system
system = ChemicalSystem(db, solution, complexation_phase_Hfo)

# Define properties
props = ChemicalProps(system)

# Specify equilibrium specs
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

# Define equilibrium solver
solver = EquilibriumSolver(specs)

# Define equilibrium conditions
conditions = EquilibriumConditions(specs)
conditions.temperature(25.0 + 273.15) # in Kelvin
conditions.pressure(1e5) # in Pa

import numpy as np
import math
pHs = np.linspace(4.0, 9.0, num=21)

metals = {"Fe": "Fe+2"}

import pandas as pd
columns = ["pH", "Hfo_sOFe+", "Hfo_wOFe+", "Hfo_wOFeOH", "Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-"]

df = pd.DataFrame(columns=columns)

# Initial Sr amount
n0 = 1e-6

def equilibrate(pH):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 2e+0, "mmol")
      state.set("Ca+2", 1e+0, "mmol")
      state.set("Fe+2", n0, "mmol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return state

print(f"  pH   mHfo_sOFe+   mHfo_wOFe+  mHfo_wOFeOH")
for pH in pHs:

      result_state = equilibrate(pH)
      mHfo_wOH    = float(result_state.speciesAmount("Hfo_wOH"))
      mHfo_wOH2   = float(result_state.speciesAmount("Hfo_wOH2+"))
      mHfo_wO     = float(result_state.speciesAmount("Hfo_wO-"))
      mHfo_sOFe   = float(result_state.speciesAmount("Hfo_sOFe+"))
      mHfo_wOFe   = float(result_state.speciesAmount("Hfo_wOFe+"))
      mHfo_wOFeOH = float(result_state.speciesAmount("Hfo_wOFeOH"))

      print(f"{pH:4.2f} {mHfo_sOFe:12.4e} {mHfo_wOFe:12.4e} {mHfo_wOFeOH:12.4e}")

      # Update dataframe with obtained values
      df.loc[len(df)] = [pH, mHfo_sOFe, mHfo_wOFe, mHfo_wOFeOH, mHfo_wOH, mHfo_wOH2, mHfo_wO]

from matplotlib import pyplot as plt
colors = ['coral', 'rosybrown', 'steelblue', 'seagreen', 'palevioletred', 'darkred', 'darkkhaki', 'cadetblue', 'indianred']

import matplotlib as mpl
from matplotlib import font_manager as fm
import os
fpath = os.path.join(mpl.get_data_path(), "texgyreadventor-regular.otf")
prop = fm.FontProperties(fname=fpath)
prop.set_size(12)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'TeX Gyre Adventor'
mpl.rcParams['font.style'] = 'normal'
mpl.rcParams['font.size'] = 14
mpl.set_loglevel("critical")

plt.figure()
ax = plt.gca()
df.plot(x="pH", y="Hfo_sOFe+", ax=ax, color=colors[0], label="Hfo_sOFe+", logy=True)
df.plot(x="pH", y="Hfo_wOFe+", ax=ax, color=colors[1], label="Hfo_wOFe+", logy=True)
df.plot(x="pH", y="Hfo_wOFeOH", ax=ax, color=colors[2], label="Hfo_wOFeOH", logy=True)
df.plot(x="pH", y="Hfo_wOH", ax=ax, color=colors[3], label="Hfo_wOH", logy=True)
df.plot(x="pH", y="Hfo_wOH2+", ax=ax, color=colors[4], label="Hfo_wOH2+", logy=True)
df.plot(x="pH", y="Hfo_wO-", ax=ax, color=colors[5], label="Hfo_wO-", logy=True)

ax.set_title("Dependence of Fe sorption on pH")
ax.set_xlabel("pH")
ax.set_ylabel("Amount of sorbed species")
ax.legend(loc="best")
ax.grid()
plt.savefig("sorbed-fe-vs-pH.png", bbox_inches='tight')
plt.close()
