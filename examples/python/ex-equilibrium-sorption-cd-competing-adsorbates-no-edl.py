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
#   • Svetlana Kyas (31 August 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O Cl Na Cd Zn Pb Cu"))
solution.setActivityModel(ActivityModelHKF())

# Define surface species lists
species_list = db.species().withAggregateState(AggregateState.Adsorbed)
species_names = extractNames(species_list)
species_names_w = [s for s in species_names if "_w" in s]
species_names_s = [s for s in species_names if "_s" in s]
species_str = ' '.join(extractNames(species_list))
species_str_w = ' '.join(species_names_w)
species_str_s = ' '.join(species_names_s)

# Create the surface
surface_Hfo = Surface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined strong site of the surface
site_Hfo_s = SurfaceSite()
site_Hfo_s.setName("Hfo_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite(site_Hfo_s)

# Defined weak site of the surface
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(species_list)

# Add specified surface as parameters for the activity model for the surface
params_site = ActivityModelSorptionParams()
params_site.surface = surface_Hfo

# Define the surface phase and set an activity model
params_site.site_tag = "_w";
hfo_w_phase = SurfacePhase(species_str_w)
hfo_w_phase.setName("Hfo_w")
hfo_w_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site))

params_site.site_tag = "_s";
hfo_s_phase = SurfacePhase(species_str_s)
hfo_s_phase.setName("Hfo_s")
hfo_s_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site))

# Create chemical system
system = ChemicalSystem(db, solution, hfo_s_phase, hfo_w_phase)

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
pHs = np.linspace(5.0, 9.0, num=17)

metals = {"Cd": "Cd+2",
          "Pb": "Pb+2",
          "Zn": "Zn+2",
          "Cu": "Cu+2"}
metal_nums = len(metals)

import pandas as pd
columns = ["pH", "Cd", "Pb", "Zn", "Cu"]
df = pd.DataFrame(columns=columns)

b_surf  = np.zeros(metal_nums)

def equilibrate(pH, metal):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 1e-1, "mol")
      state.set("Na+" , 1e-1, "mol")
      state.set("Cd+2", 5e-7, "mol")
      if metal == "Pb":
            state.set("Pb+2", 5e-5, "mol")
      elif metal == "Cu":
            state.set("Cu+2", 5e-5, "mol")
      elif metal == "Zn":
            state.set("Zn+2", 5e-5, "mol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return res, state

print(f"  pH    % sorbed Cd (only Cd)  % sorbed Cd (Cd + Pb)   % sorbed Cd (Cd + Zn)   % sorbed Cd (Cd + Cu)")

for pH in pHs:

      print(f"{pH:4.2f}", end="\t")
      i = 0
      for metal in metals:

            res, state = equilibrate(pH, metal)

            # If the equilibrium calculations didn't succeed, continue to the next condition
            if res.optima.succeeded:
                  # Update properties
                  props.update(state)

                  # Fetch amount of sorbed and dissolved mineral
                  b_total = float(props.elementAmount("Cd"))
                  b_surf[i]  = float(props.elementAmountInPhase("Cd", "Hfo_w") + props.elementAmountInPhase("Cd", "Hfo_s")) / b_total * 100
                  print(f"{b_surf[i]:20.4f}", end="\t")
            i += 1
      print("")

      # Update dataframe with obtained values
      df.loc[len(df)] = [pH] + list(b_surf)

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
df.plot(x="pH", y="Cd", ax=ax, color=colors[1], label="5e-7 mol Cd")
df.plot(x="pH", y="Cu", ax=ax, color=colors[0], marker="v", label="5e-7 mol Cd + 5e-5 mol Cu")
df.plot(x="pH", y="Zn", ax=ax, color=colors[2], marker="o", label="5e-7 mol Cd + 5e-5 mol Zn")
df.plot(x="pH", y="Pb", ax=ax, color=colors[3], marker="s", label="5e-7 mol Cd + 5e-5 mol Pb")
ax.set_title("Cd sorption competition with Pb, Zn, and Cu vs pH (no EDL)")
ax.set_xlabel("pH")
ax.set_ylabel("% of sorbed metal")
ax.legend(loc="best")
ax.grid()
plt.savefig("competing-adsorbates-vs-pH-no-edl.png", bbox_inches='tight')
plt.close()
