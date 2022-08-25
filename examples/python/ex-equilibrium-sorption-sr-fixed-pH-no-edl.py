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
solution = AqueousPhase(speciate("H O Cl Ca Sr"))
solution.setActivityModel(ActivityModelHKF())

# Define surface species list
list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH"
all_species = db.species().withAggregateState(AggregateState.Adsorbed)
list = all_species.withNames(list_str)

list_str_s = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2"
list_str_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH"

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
conditions.pressure(1.0, "bar")
conditions.pH(6.0)

# Define initial equilibrium state
state = ChemicalState(system)
state.set("H2O" , 1.00, "kg")
state.set("Cl-" , 2e+0, "mmol")
state.set("Ca+2", 1e+0, "mmol")
state.set("Sr+2", 1e-6, "mmol")
state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

# Equilibrate given initial state with input conditions
res = solver.solve(state, conditions)

print("*******************************************\n"
      "After equilibration:\n"
      "*******************************************")
print("succeed       = ", res.optima.succeeded)
print("solutionstate = \n", state)

# Evaluate aqueous properties
aqprops = AqueousProps(state)
print("I  = %e mol/kgw" % float(aqprops.ionicStrength()))
print("pH = %f" % float(aqprops.pH()))

# Evaluate surface and sites properties
sprops = SurfaceProps(surface_Hfo, state)
print(sprops)