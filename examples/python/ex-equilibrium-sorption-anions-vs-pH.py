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
#   • Svetlana Kyas (30 July 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O S C P"))
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
pHs = np.linspace(2.0, 9.0, num=29)

# anions = {"SO4-2": ["S", ["Hfo_wOHSO4-2", "Hfo_wSO4-"]],
#           "PO4-3": ["P", ["Hfo_wH2PO4", "Hfo_wHPO4-", "Hfo_wPO4-2"]],
#           "CO3-2": ["C", ["Hfo_wCO3-", "Hfo_wHCO3"]]}
anions = {"SO4-2": ["S", ["Hfo_wOHSO4-2", "Hfo_wSO4-"]]}

import pandas as pd
columns = ["Anion", "pH", "%"]
df = pd.DataFrame(columns=columns)

# Initial Sr amount
n0 = 1e-3

def equilibrate(pH, anion):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      # state.set("Cl-" , 2e+0, "mmol")
      # state.set("Ca+2", 1e+0, "mmol")
      state.set(anion, n0, "mol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      # If the equilibrium calculations didn't succeed, continue to the next condition
      if res.optima.succeeded:
            # Update properties
            props.update(state)

            # Fetch amount of sorbed and dissolved mineral
            n_aq = float(props.elementAmountInPhase(anions[anion][0], "AqueousPhase"))
            n_surf = 0
            for s in anions[anion][1]:
                  #print("surf species", s)
                  n_surf_tmp = float(state.speciesAmount(s))
                  n_surf += n_surf_tmp
            # print("n_surf", n_surf)
            # print("n_aq", n_aq)
            # print("n_total", n_aq + n_surf)
            # input()

            return n_surf, n_aq
      else:
            return math.nan, math.nan

for anion in anions:
      print(f"  pH   % sorbed {anion}   % dissolved {anion}")

      for pH in pHs:

            result = equilibrate(pH, anion)
            n_surf = result[0] / n0 * 100
            n_aq = result[1] / n0 * 100
            print(f"{pH:4.2f} {n_surf:12.4f} {n_aq:16.4f}")

            # Update dataframe with obtained values
            df.loc[len(df)] = [anion, pH, n_surf]

from matplotlib import pyplot as plt
colors = ['coral', 'rosybrown', 'steelblue', 'seagreen', 'palevioletred', 'darkred', 'darkkhaki', 'cadetblue', 'indianred']

plt.figure()
df_anion = df[df["Anion"] == list(anions)[0]] # fetch the columns with Pb+2
ax = df_anion.plot(x="pH", y="%", color=colors[0], label=list(anions.keys())[0])
ax.set_title("Dependence of anion sorption on pH")
ax.set_xlabel("pH")
ax.set_ylabel("% of sorbed anion")
for idx, anion in enumerate(anions):
      if idx:
            df_anion = df[df["Anion"] == anion] # fetch the columns with other anions
            df_anion.plot(x="pH", y="%", ax=ax, color=colors[idx+1], label=anion)
ax.legend(loc="best")
ax.grid()
plt.savefig("sorbed-anions-vs-pH.png", bbox_inches='tight')
plt.close()
