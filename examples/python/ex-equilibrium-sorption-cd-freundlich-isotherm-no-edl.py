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
#   • Svetlana Kyas (9 June 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase(speciate("H O Cl Ca Cd"))
solution.setActivityModel(ActivityModelHKF())

# Define ion exchange species list
list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ " \
           "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+"
all_species = db.species().withAggregateState(AggregateState.Adsorbed)
list = all_species.withNames(list_str)

list_str_s = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+"
list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+"

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined sites of the complexation surface
surface_Hfo.addSite("Hfo_s", "_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(list)
print(surface_Hfo)

# Add specified surface as parameters for the activity model for the complexation surface
params_site = ActivityModelSurfaceComplexationSiteParams()
params_site.surface = surface_Hfo

# Define surface complexation phases and set an activity model
params_site.site_tag = "_w";
hfo_w_phase = SurfaceComplexationPhase(list_str_w)
hfo_w_phase.setName("Hfo_w")
hfo_w_phase.setActivityModel(ActivityModelSurfaceComplexationSiteNoDDL(params_site))
params_site.site_tag = "_s";
hfo_s_phase = SurfaceComplexationPhase(list_str_s)
hfo_s_phase.setName("Hfo_s")
hfo_s_phase.setActivityModel(ActivityModelSurfaceComplexationSiteNoDDL(params_site))

# Create chemical system
system = ChemicalSystem(db, solution, hfo_w_phase, hfo_s_phase)

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
conditions.temperature(25.0, "celsius")
conditions.pressure(1.0, "atm")
conditions.pH(6.0)

import numpy as np

# Mass of given, sorbed and dissolved Cd (in mg)
mCds = np.linspace(1, 40, num=40)
mCdsurf = []
mCdaq = []

def equilibrate(mCd):

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 2e+0, "mmol")
      state.set("Ca+2", 1e+0, "mmol")
      state.set("Cd+2",  mCd, "mg")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return state

print(f" Input m(Cd), mg  Sorbed m(Cd), mg  Dissolved m(Cd), mg Kd=sI/cI [-]")
for mCd in mCds:

      state = equilibrate(mCd)

      # Update properties
      props.update(state)

      # Fetch total, dissolved, and sorbed mass of Cd (converted in mg)
      massCd_aq   = 1e6*float(props.elementMassInPhase("Cd", "AqueousPhase"))
      massCd_surf = 1e6*float(props.elementMassInPhase("Cd", "Hfo_w") + props.elementMassInPhase("Cd", "Hfo_s"))

      mCdaq.append(massCd_aq)
      mCdsurf.append(massCd_surf)

      print(f"{mCd:16.2f} {massCd_surf:17.4e} {massCd_aq:20.4e} {massCd_surf/massCd_aq:12.4f}")

from matplotlib import pyplot as plt
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

plt.plot(mCdaq, mCdsurf)
plt.title(r"Freundlich sorption isotherm (no EDL)")
plt.xlabel('Dissolved mass of Cd [mg]')
plt.ylabel('Sorbed mass of Cd [mg]')
plt.savefig("sorbed-cd-mass-vs-dissolved-cd-mass-no-edl.png", bbox_inches='tight')
plt.close()

plt.plot(mCdaq, np.array(mCdsurf) / np.array(mCdaq))
plt.title("Distribution coefficient \n of Freundlich sorption isotherm (no EDL)")
plt.xlabel('Dissolved mass of Cd [mg]')
plt.ylabel(r'K$_d$ = s$_I$ / c$_I$ [-]')
plt.savefig("distribution-coeff-vs-dissolved-cd-mass-no-edl.png", bbox_inches='tight')
plt.close()
