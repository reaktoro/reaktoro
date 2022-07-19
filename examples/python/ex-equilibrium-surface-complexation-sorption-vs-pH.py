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
solution = AqueousPhase(speciate("H O Cl Ca Sr Cd Zn Pb Cu Fe Ba"))
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
conditions.temperature(25.0, "celsius")
conditions.pressure(1.0, "atm")

import numpy as np
pHs = np.linspace(2.0, 9.0, num=15)
b_surf = []
b_aq = []

metal = {"Pb": "Pb+2",
         "Cd": "Cd+2",
         "Zn": "Zn+2",
         "Cu": "Cu+2",
         "Sr": "Sr+2",
         "Fe": "Fe+3",
         "Ba": "Ba+2",
         "Ca": "Ca+2"}

#
#metal_b = "Pb"
#metal_b = "Cd"
#metal_b = "Zn"
#metal_b = "Cu"
#metal_b = "Sr"
#metal_b = "Fe"
#metal_b = "Ba"
metal_b = "Ca"

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
      state.set(metal[metal_b], n0, "mmol")
      state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
      state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

      # Equilibrate given initial state with input conditions
      res = solver.solve(state, conditions)

      # Update properties
      props.update(state)

      # Fetch amount of sorbed and dissolved mineral
      b_aq    = float(props.elementAmountInPhase(metal_b, "AqueousPhase"))
      b_surf  = float(props.elementAmountInPhase(metal_b, "SurfaceComplexationPhase"))

      return b_surf, b_aq

print(f"        pH  % sorbed {metal_b}")
for pH in pHs:

      result = equilibrate(pH)
      b_total = float(props.elementAmount(metal_b))
      b_surf.append(result[0] / b_total * 100)
      b_aq.append(result[1] / b_total * 100)

      print(f"{pH:6.4e} {b_surf[-1]:12.4f}")

from matplotlib import pyplot as plt
plt.plot(pHs, b_surf)
plt.title(f"Dependence of {metal_b} sorption on pH")
plt.xlabel('pH [-]')
plt.ylabel(f'Sorbed {metal_b} [%]')
plt.grid()
plt.savefig("sorbed-" + metal_b + "-vs-pH.png", bbox_inches='tight')
plt.close()
