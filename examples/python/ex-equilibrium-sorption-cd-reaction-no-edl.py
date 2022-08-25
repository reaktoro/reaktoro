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
solution = AqueousPhase(speciate("H O Cl Ca Cd"))
solution.setActivityModel(ActivityModelHKF())

# Define ion exchange species list
list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ " \
           "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+"
all_species = db.species().withAggregateState(AggregateState.Adsorbed)
list = all_species.withNames(list_str)

list_str_s = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+"
list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+"

# Create the surface
surface_Hfo = Surface("Hfo")
surface_Hfo.setSpecificSurfaceArea(60, "m2/g").setMass(4.45, "g")

# Defined sites of the surface
surface_Hfo.addSite("Hfo_s", "_s").setAmount(0.025e-3, "mol")
surface_Hfo.addSite("Hfo_w", "_w").setAmount(1e-3, "mol")

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(list)
print(surface_Hfo)

# Add specified surface as parameters for the activity model for the surface
params_site = ActivityModelSorptionParams()
params_site.surface = surface_Hfo

# Define surface phases and set an activity model
params_site.site_tag = "_w";
hfo_w_phase = SurfacePhase(list_str_w)
hfo_w_phase.setName("Hfo_w")
hfo_w_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site))
params_site.site_tag = "_s";
hfo_s_phase = SurfacePhase(list_str_s)
hfo_s_phase.setName("Hfo_s")
hfo_s_phase.setActivityModel(ActivityModelSorptionNoDDL(params_site))

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
mCd = 0
delta_mCd = 1e-1
num_steps = 20
mCdsurf = []
mCdaq = []
mCdtot = []

# Define initial equilibrium state
state = ChemicalState(system)
state.set("H2O" , 1.00, "kg")
state.set("Cl-" , 2e+0, "mmol")
state.set("Ca+2", 1e+0, "mmol")
state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

def equilibrate():

      state.add("Cd+2", delta_mCd, "mmol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      return state

print(f"Input m(Cd), mmol   Sorbed m(Cd), mmol   Dissolved m(Cd), mmol   Kd=sI/cI [-]    Kd = Sorbed Cd / Total Cd [-]")
for step in range(0, num_steps):

      mCd += delta_mCd
      state = equilibrate()

      # Update properties
      props.update(state)

      # Fetch total, dissolved, and sorbed mass of Cd
      massCd_aq   = 1e3*float(props.elementAmountInPhase("Cd", "AqueousPhase"))
      massCd_surf = 1e3*float(props.elementMassInPhase("Cd", "Hfo_w") + props.elementMassInPhase("Cd", "Hfo_s"))
      massCd_tot  = 1e3*float(props.elementAmount("Cd"))

      mCdaq.append(massCd_aq)
      mCdsurf.append(massCd_surf)
      mCdtot.append(massCd_tot)

      print(f"{mCd:17.2e} {mCdsurf[-1]:20.4e} {mCdaq[-1]:23.4e} {mCdsurf[-1]/mCdaq[-1]:14.4e} {mCdsurf[-1]/mCdtot[-1]:32.4e} ")

from matplotlib import pyplot as plt
plt.plot(mCdaq, mCdsurf, color='coral')
plt.title(r"Freundlich sorption isotherm (no EDL)")
plt.xlabel('Dissolved mass of Cd [mmol]')
plt.ylabel('Sorbed mass of Cd [mmol]')
plt.savefig("sorbed-cd-mass-vs-cd-reaction-no-edl.png", bbox_inches='tight')
plt.close()

plt.plot(mCdaq, np.array(mCdsurf) / np.array(mCdaq), color='coral')
plt.title("Distribution coefficient \n of Freundlich sorption isotherm \n (no EDL) ")
plt.xlabel('Dissolved mass of Cd [mmol]')
plt.ylabel(r'K$_d$ = s$_I$ / c$_I$ [-]')
plt.savefig("kd-vs-cd-reaction-no-edl.png", bbox_inches='tight')
plt.close()
