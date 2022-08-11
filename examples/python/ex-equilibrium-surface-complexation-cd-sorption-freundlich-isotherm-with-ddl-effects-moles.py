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
list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOCd+ Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOCd+"
all_species = db.species().withAggregateState(AggregateState.Adsorbed)
list = all_species.withNames(list_str)

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
surface_Hfo.addSurfaceSpecies(list)

# Add specified surface as parameters for the activity model for the complexation surface
params = ActivityModelSurfaceComplexationParams()
params.surface = surface_Hfo

# Define surface complexation phase and set an activity model
complexation_phase_Hfo = SurfaceComplexationPhase(list_str)
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params))

# Create chemical system
system = ChemicalSystem(db, solution, complexation_phase_Hfo)

props = ChemicalProps(system)
aqprops = AqueousProps(system)
surfprops = ComplexationSurfaceProps(surface_Hfo, system)

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
mCds = np.linspace(10, 100, num=21)
mCdsurf = []
mCdaq = []
mCdtot = []


def equilibrate(mCd):

      # Define initial equilibrium state
      state = ChemicalState(system)
      state.set("H2O" , 1.00, "kg")
      state.set("Cl-" , 2e+0, "mmol")
      state.set("Ca+2", 1e+0, "mmol")
      state.set("Cd+2",  mCd, "mmol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # print(state.speciesAmount("Cd+2"))
      # input()
      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      # Update properties
      props.update(state)

      # Fetch total, dissolved, and sorbed mass of Cd
      massCd_aq   = float(props.elementAmountInPhase("Cd", "AqueousPhase"))
      massCd_surf = float(props.elementAmountInPhase("Cd", "SurfaceComplexationPhase"))
      massCd_tot  = float(props.elementAmount("Cd"))

      return massCd_aq, massCd_surf, massCd_tot

print(f"Input(Cd) [mol]   Total(Cd) [mol]   Sorbed(Cd) [mol]   Dissolved(Cd) [mol]   Sorbed(Cd)/Dissolved(Cd) [-]   Sorbed(Cd)/Total(Cd) [-]")
for mCd in mCds:

      result = equilibrate(mCd)
      mCdaq.append(result[0])
      mCdsurf.append(result[1])
      mCdtot.append(result[2])

      print(f"{mCd*1e-3:15.4e} {mCdtot[-1]:17.4e} {mCdsurf[-1]:18.4e} {mCdaq[-1]:21.4e} {mCdsurf[-1]/mCdaq[-1]:30.4e} {mCdsurf[-1]/mCdtot[-1]:26.4e} ")

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

plt.plot(mCdaq, mCdsurf, color='coral')
plt.title(r"Freundlich sorption isotherm (with electrostatic effects)")
plt.xlabel('Dissolved mass of Cd [mol]')
plt.ylabel('Sorbed mass of Cd [mol]')
plt.grid()
plt.savefig("sorbed-cd-mass-vs-dissolved-cd-mass-with-ddl-effects-moles.png", bbox_inches='tight')
plt.close()

plt.plot(mCdaq, np.array(mCdsurf) / np.array(mCdaq), color='coral')
plt.title("Distribution coefficient of Freundlich sorption isotherm \n (with electrostatic effects)")
plt.xlabel('Dissolved mass of Cd [mol]')
plt.ylabel(r'K$_d$ = s$_I$ / c$_I$ [-]')
plt.grid()
plt.savefig("distribution-coeff-vs-dissolved-cd-mass-with-ddl-effects-moles.png", bbox_inches='tight')
plt.close()
