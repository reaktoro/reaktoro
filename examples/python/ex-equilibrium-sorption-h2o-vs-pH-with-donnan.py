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
solution = AqueousPhase(speciate("H O Na Cl"))
solution.setActivityModel(ActivityModelHKF())

# Define surface complexation  species list
list_str = "Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_wOH Hfo_wOH2+ Hfo_wO-";
slist = db.species().withAggregateState(AggregateState.Adsorbed)
list = slist.withNames(StringList(list_str))

list_str_s = "Hfo_sOH Hfo_sOH2+ Hfo_sO-"
list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO-"

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined sites of the complexation surface
surface_Hfo.addSite("Hfo_s", "_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(list)

# Add specified surface as parameters for the activity model for the complexation surface
params_site = ActivityModelSurfaceComplexationSiteParams()
params_site.surface = surface_Hfo

# Define surface complexation phase and set an activity model
params_site.site_tag = "_w";
complexation_phase_Hfo_w = SurfaceComplexationPhase(list_str_w)
complexation_phase_Hfo_w.setName("Hfo_w")
complexation_phase_Hfo_w.setActivityModel(ActivityModelSurfaceComplexationSiteWithDDL(params_site))

params_site.site_tag = "_s";
complexation_phase_Hfo_s = SurfaceComplexationPhase(list_str_s)
complexation_phase_Hfo_s.setName("Hfo_s")
complexation_phase_Hfo_s.setActivityModel(ActivityModelSurfaceComplexationSiteWithDDL(params_site))

# Define the DDL phase
ddl_Hfo = DoubleLayerPhase(speciate("H O Na Cl"))
ddl_Hfo.named("DoubleLayerPhase")
params_dll = ActivityModelDDLParams()
ddl_Hfo.setActivityModel(chain(ActivityModelHKF(), ActivityModelDonnanDDL(ActivityModelDDLParams())))

# Create chemical system
system = ChemicalSystem(db, solution, complexation_phase_Hfo_w, complexation_phase_Hfo_s, ddl_Hfo)

# Define properties
props = ChemicalProps(system)
aprops = AqueousProps(system)
site_w_props = ComplexationSurfaceSiteProps(surface_Hfo.sites()["_w"], system)
site_s_props = ComplexationSurfaceSiteProps(surface_Hfo.sites()["_s"], system)
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
pHs = np.linspace(4.0, 10.0, num=31)

import pandas as pd
columns = ["pH",
           "Hfo_sOH", "Hfo_sOH2+", "Hfo_sO-", "Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-",
           "xHfo_sOH", "xHfo_sOH2+", "xHfo_sO-", "xHfo_wOH", "xHfo_wOH2+", "xHfo_wO-"]

df = pd.DataFrame(columns=columns)

# Initial Sr amount
nNaCl = 0.0
def equilibrate(pH):

      # Set pH
      conditions.pH(pH)

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O", 1.00, "kg")
      state.set("Cl-", nNaCl, "mol")
      state.set("Na+", nNaCl, "mol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return state, res

print(f"  pH   x(Hfo_sOH) x(Hfo_sOH2+)   x(Hfo_sO-)   x(Hfo_wOH) x(Hfo_wOH2+)   x(Hfo_wO-) succeed I")
for pH in pHs:

      result_state, res = equilibrate(pH)

      mHfo_wOH  = float(result_state.speciesAmount("Hfo_wOH"))
      mHfo_wOH2 = float(result_state.speciesAmount("Hfo_wOH2+"))
      mHfo_wO   = float(result_state.speciesAmount("Hfo_wO-"))
      mHfo_sOH  = float(result_state.speciesAmount("Hfo_sOH"))
      mHfo_sOH2 = float(result_state.speciesAmount("Hfo_sOH2+"))
      mHfo_sO   = float(result_state.speciesAmount("Hfo_sO-"))

      site_w_props.update(result_state)
      site_s_props.update(result_state)
      aprops.update(result_state)

      xHfo_wOH  = float(site_w_props.speciesFraction("Hfo_wOH"))
      xHfo_wOH2 = float(site_w_props.speciesFraction("Hfo_wOH2+"))
      xHfo_wO   = float(site_w_props.speciesFraction("Hfo_wO-"))
      xHfo_sOH  = float(site_s_props.speciesFraction("Hfo_sOH"))
      xHfo_sOH2 = float(site_s_props.speciesFraction("Hfo_sOH2+"))
      xHfo_sO   = float(site_s_props.speciesFraction("Hfo_sO-"))

      print(f"{pH:4.1f} {xHfo_sOH:12.4e} {xHfo_sOH2:12.4e} {xHfo_sO:12.4e} {xHfo_wOH:12.4e} {xHfo_wOH2:12.4e} {xHfo_wO:12.4e} {res.optima.succeeded:6d} {float(aprops.ionicStrength()):6.2f}")

      # Update dataframe with obtained values
      df.loc[len(df)] = [pH,
                         mHfo_sOH, mHfo_sOH2, mHfo_sO, mHfo_wOH, mHfo_wOH2, mHfo_wO,
                         xHfo_sOH, xHfo_sOH2, xHfo_sO, xHfo_wOH, xHfo_wOH2, xHfo_wO]

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
df.plot(x="pH", y="Hfo_sOH", ax=ax, color=colors[0], label=r"Hfo_sOH", logy=True)
df.plot(x="pH", y="Hfo_sOH2+", ax=ax, color=colors[1], label=r"Hfo_sOH$_2^+$", logy=True)
df.plot(x="pH", y="Hfo_sO-", ax=ax, color=colors[2], label=r"Hfo_sO_$-$", logy=True)
df.plot(x="pH", y="Hfo_wOH", ax=ax, color=colors[3], label=r"Hfo_wOH", logy=True)
df.plot(x="pH", y="Hfo_wOH2+", ax=ax, color=colors[4], label=r"Hfo_wOH$_2^+$", logy=True)
df.plot(x="pH", y="Hfo_wO-", ax=ax, color=colors[5], label=r"Hfo_wO$^-$", logy=True)

ax.set_title("Dependence of water sorption on pH")
ax.set_xlabel("pH")
ax.set_ylabel("Amount of sorbed species")
ax.legend(loc="best")
ax.grid()
plt.savefig(f"sorbed-amounts-h2o-vs-pH-with-esurf-NaCl-{nNaCl}.png", bbox_inches='tight')
plt.close()

plt.figure()
ax = plt.gca()
df.plot(x="pH", y="xHfo_sOH"  , ax=ax, color=colors[0], label=r"Hfo_sOH"      )
df.plot(x="pH", y="xHfo_sOH2+", ax=ax, color=colors[1], label=r"Hfo_sOH$_2^+$"  )
df.plot(x="pH", y="xHfo_sO-"  , ax=ax, color=colors[2], label=r"Hfo_sO$^-$")
df.plot(x="pH", y="xHfo_wOH"  , ax=ax, color=colors[3], label=r"Hfo_wOH"   )
df.plot(x="pH", y="xHfo_wOH2+", ax=ax, color=colors[4], label=r"Hfo_wOH$_2^+$")
df.plot(x="pH", y="xHfo_wO-"  , ax=ax, color=colors[5], label=r"Hfo_wO$^-$"   )

ax.set_title("Dependence of water sorption on pH")
ax.set_xlabel("pH")
ax.set_ylabel("Fraction of sorbed species")
ax.legend(loc="best")
ax.grid()
plt.savefig(f"sorbed-fractions-h2o-vs-with-esurf-pH-NaCl-{nNaCl}.png", bbox_inches='tight')
plt.close()