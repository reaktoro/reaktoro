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

# Define surface complexation species lists
species_list = db.species().withAggregateState(AggregateState.Adsorbed)
species_names = extractNames(species_list)
species_names_w = [s for s in species_names if "_w" in s]
species_names_s = [s for s in species_names if "_s" in s]
species_str = ' '.join(extractNames(species_list))
species_str_w = ' '.join(species_names_w)
species_str_s = ' '.join(species_names_s)

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined strong and weak sites of the complexation surface
surface_Hfo.addSite("Hfo_s", "_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(species_list)

# Add specified surface as parameters for the activity model for the complexation surface
params_site = ActivityModelSurfaceComplexationSiteParams()
params_site.surface = surface_Hfo

# Define surface complexation phase and set an activity model
params_site.site_tag = "_w";
hfo_w_phase = SurfaceComplexationPhase(species_str_w)
hfo_w_phase.setName("Hfo_w")
hfo_w_phase.setActivityModel(ActivityModelSurfaceComplexationSiteNoDDL(params_site))

params_site.site_tag = "_s";
hfo_s_phase = SurfaceComplexationPhase(species_str_s)
hfo_s_phase.setName("Hfo_s")

# Create chemical system
system = ChemicalSystem(db, solution, hfo_s_phase, hfo_w_phase)

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
pHs = np.linspace(4.0, 9.0, num=21)

metals = {"Fe": "Fe+2"}

import pandas as pd
columns = ["pH", "Hfo_wOFe+", "Hfo_wOFeOH", "Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-",
                 "xHfo_wOFe+", "xHfo_wOFeOH", "xHfo_wOH", "xHfo_wOH2+", "xHfo_wO-"]

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

print(f"  pH     xHfo_wOH   xHfo_wOH2+     xHfo_wO-   xHfo_wOFe+  xHfo_wOFeOH            I            Z        sigma")
for pH in pHs:

      result_state = equilibrate(pH)

      mHfo_wOH    = float(result_state.speciesAmount("Hfo_wOH"))
      mHfo_wOH2   = float(result_state.speciesAmount("Hfo_wOH2+"))
      mHfo_wO     = float(result_state.speciesAmount("Hfo_wO-"))
      mHfo_wOFe   = float(result_state.speciesAmount("Hfo_wOFe+"))
      mHfo_wOFeOH = float(result_state.speciesAmount("Hfo_wOFeOH"))

      site_w_props.update(result_state)
      site_s_props.update(result_state)
      aprops.update(result_state)

      xHfo_wOH    = float(site_w_props.speciesFraction("Hfo_wOH"))
      xHfo_wOH2   = float(site_w_props.speciesFraction("Hfo_wOH2+"))
      xHfo_wO     = float(site_w_props.speciesFraction("Hfo_wO-"))
      xHfo_wOFe   = float(site_w_props.speciesFraction("Hfo_wOFe+"))
      xHfo_wOFeOH = float(site_w_props.speciesFraction("Hfo_wOFeOH"))

      Z_s = site_s_props.charge();
      Z_w = site_w_props.charge();
      Z = float(Z_s + Z_w)
      sigma = float(site_w_props.sigma(Z_w) + site_s_props.sigma(Z_s))

      print(f"{pH:4.1f} {xHfo_wOH:12.4e} {xHfo_wOH2:12.4e} {xHfo_wO:12.4e} {xHfo_wOFe:12.4e} {xHfo_wOFeOH:12.4e} {float(aprops.ionicStrength()):12.4e} {Z:12.4e} {sigma:12.4e}")

      # Update dataframe with obtained values
      df.loc[len(df)] = [pH, mHfo_wOFe, mHfo_wOFeOH, mHfo_wOH, mHfo_wOH2, mHfo_wO,
                             xHfo_wOFe, xHfo_wOFeOH, xHfo_wOH, xHfo_wOH2, xHfo_wO]

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
df.plot(x="pH", y="Hfo_wOFe+",  ax=ax, color=colors[1], label=r"Hfo_wOFe$^+$",  logy=True)
df.plot(x="pH", y="Hfo_wOFeOH", ax=ax, color=colors[2], label="Hfo_wOFeOH",     logy=True)
df.plot(x="pH", y="Hfo_wOH",    ax=ax, color=colors[3], label="Hfo_wOH",        logy=True)
df.plot(x="pH", y="Hfo_wOH2+",  ax=ax, color=colors[4], label=r"Hfo_wOH$_2^+$", logy=True)
df.plot(x="pH", y="Hfo_wO-",    ax=ax, color=colors[5], label=r"Hfo_wO$^-$",    logy=True)

ax.set_title("Dependence of Fe sorption on pH (noo EDL)")
ax.set_xlabel("pH")
ax.set_ylabel("Amounts of sorbed species")
ax.legend(loc="best")
ax.grid()
plt.savefig("sorbed-amounts-fe-h2o-vs-pH-no-edl.png", bbox_inches='tight')
plt.close()

plt.figure()
ax = plt.gca()
df.plot(x="pH", y="xHfo_wOFe+",  ax=ax, color=colors[1], label=r"Hfo_wOFe$^+$")
df.plot(x="pH", y="xHfo_wOFeOH", ax=ax, color=colors[2], label="Hfo_wOFeOH")
df.plot(x="pH", y="xHfo_wOH",    ax=ax, color=colors[3], label="Hfo_wOH")
df.plot(x="pH", y="xHfo_wOH2+",  ax=ax, color=colors[4], label=r"Hfo_wOH$_2^+$")
df.plot(x="pH", y="xHfo_wO-",    ax=ax, color=colors[5], label=r"Hfo_wO$^-$")

ax.set_title("Dependence of Fe sorption on pH (no EDL)")
ax.set_xlabel("pH")
ax.set_ylabel("Fractions of sorbed species")
ax.legend(loc="best")
ax.grid()
plt.savefig("sorbed-fractions-fe-h2o-vs-pH-no-edl.png", bbox_inches='tight')
plt.close()
