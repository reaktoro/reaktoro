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
solution = AqueousPhase(speciate("H O Cl S Cd P Na C F"))
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
pHs = np.linspace(6.0, 8.0, num=13)

solutions = {"Cd": "Cd+2",
             "SO4": "SO4-",
             "Cl": "Cl-",
             "PO4": "PO4-3",
             "CO3": "CO3-2",
             "F": "F-"}
solution_nums = len(solutions)

import pandas as pd
columns = ["pH", "Cd", "SO4", "Cl", "PO4", "CO3", "F"]
df = pd.DataFrame(columns=columns)

b_surf  = np.zeros(solution_nums)

def equilibrate(pH, anion):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 7e-1, "mol")
      state.set("Na+" , 7e-1, "mol")
      state.set("Cd+2", 5e-7, "mol")
      if anion == "SO4":
            state.set("SO4-2", 2e-1, "mol")
      elif anion == "Cl":
            state.set("Cl-", 5e-1, "mol")
      elif anion == "PO4":
            state.set("PO4-3", 7e-1, "mol")
      elif anion == "CO3":
            state.set("CO3-2", 1e-1, "mol")
      elif anion == "F":
            state.set("F-", 9e-1, "mol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return res, state

print(f"  pH    % sorbed Cd (no anions)   % sorbed Cd (with SO4-)  % sorbed Cd (with Cl-)   % sorbed Cd (with PO4-3)  % sorbed Cd (with CO3-2) % sorbed Cd (with F-)")

for pH in pHs:

      print(f"{pH:4.2f}", end="\t")
      i = 0
      for anion in solutions:

            res, state = equilibrate(pH, anion)

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
df.plot(x="pH", y="Cd", ax=ax, color=colors[5], label="5e-7 moles Cd", marker="8")
df.plot(x="pH", y="SO4", ax=ax, color=colors[0], label="5e-7 moles Cd + 0.2 mole SO4-2", marker="^")
df.plot(x="pH", y="Cl", ax=ax, color=colors[1], label="5e-7 moles Cd + 0.5 mole Cl-", marker="o")
#df.plot(x="pH", y="PO4", ax=ax, color=colors[2], label="0.7 mole PO4-3", marker="+")
#df.plot(x="pH", y="CO3", ax=ax, color=colors[3], label="0.1 mole CO3-2", marker="*")
#df.plot(x="pH", y="F", ax=ax, color=colors[4], label="0.9 mole F-", marker="s")
ax.set_title("Effect of anions on Cd sorption vs  pH (no EDL)")
ax.set_xlabel("pH")
ax.set_ylabel("% of sorbed Cd")
ax.legend(loc="best")
ax.grid()
plt.savefig("cd-sorption-solution-effect-vs-pH-no-edl.png", bbox_inches='tight')
plt.close()
