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
#   • Svetlana Kyas (9 August 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O Cl Ca Sr Cd Zn Pb Cu Fe"))
solution.setActivityModel(ActivityModelHKF())

# Define surface complexation  species list
species_list = db.species().withAggregateState(AggregateState.Adsorbed)
species_str = ' '.join(extractNames(species_list))

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined weak site of the complexation surface
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Defined strong site of the complexation surface
site_Hfo_s = ComplexationSurfaceSite()
site_Hfo_s.setName("Hfo_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite(site_Hfo_s)

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(species_list)

# Add specified surface as parameters for the activity model for the complexation surface
params = ActivityModelSurfaceComplexationParams()
params.surface = surface_Hfo

# Define surface complexation phase and set an activity model
complexation_phase_Hfo = SurfaceComplexationPhase(species_str)
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params))

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
pHs = np.linspace(2.0, 9.0, num=29)

metals = {"Pb": "Pb+2",
          "Cd": "Cd+2",
          "Zn": "Zn+2",
          "Cu": "Cu+2",
          "Sr": "Sr+2",
          "Fe": "Fe+2"}

import pandas as pd
columns = ["Metal", "pH", "%"]
df = pd.DataFrame(columns=columns)

# Initial Sr amount
n0 = 1e-6

def equilibrate(pH, metal):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 2e+0, "mmol")
      state.set("Ca+2", 1e+0, "mmol")
      state.set(metals[metal], n0, "mmol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      # If the equilibrium calculations didn't succeed, continue to the next condition
      if res.optima.succeeded:
            # Update properties
            props.update(state)

            # Fetch amount of sorbed and dissolved mineral
            b_aq    = float(props.elementAmountInPhase(metal, "AqueousPhase"))
            b_surf  = float(props.elementAmountInPhase(metal, "SurfaceComplexationPhase"))

            return b_surf, b_aq
      else:
            return math.nan, math.nan

for metal in metals:
      print(f"pH    % sorbed {metal}   % dissolved {metal}")

      for pH in pHs:

            result = equilibrate(pH, metal)
            b_total = float(props.elementAmount(metal))
            if b_total != 0:
                  b_surf = result[0] / b_total * 100
                  b_aq = result[1] / b_total * 100
                  print(f"{pH:4.2f} {b_surf:12.4f} {b_aq:16.4f}")
            else:
                  b_surf, b_aq = result[0], result[1]

            # Update dataframe with obtained values
            df.loc[len(df)] = [metals[metal], pH, b_surf]

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
df_metal = df[df["Metal"] == list(metals.values())[0]] # fetch the columns with Pb+2
ax = df_metal.plot(x="pH", y="%", color=colors[0], label=list(metals.keys())[0])
ax.set_xlabel("pH")
ax.set_ylabel("% of sorbed metal")
ax.set_title("Dependence of the metal sorption on pH \n (with electrostatic effects)")
for idx, metal in enumerate(metals):
      if idx:
            df_metal = df[df["Metal"] == metals[metal]] # fetch the columns with other metals
            df_metal.plot(x="pH", y="%", ax=ax, color=colors[idx+1], label=metal)
ax.legend(loc="best")
ax.grid()
plt.savefig("sorbed-metals-vs-pH-ddl-effects.png", bbox_inches='tight')
plt.close()
