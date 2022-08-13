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
#   • Svetlana Kyas (8 June 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------


from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2 Mg+2 Sr+2 Cd+2")
solution.setActivityModel(ActivityModelHKF())

# Define ion exchange species list
list_str = "Hfo_sOH Hfo_sOCd+ Hfo_sOHSr+2 Hfo_wOH Hfo_wOCd+ Hfo_wOSr+ Hfo_wOSrOH"
all_species = db.species().withAggregateState(AggregateState.Adsorbed)
list = all_species.withNames(list_str)

# Create complexation surface
surface_Hfo = ComplexationSurface("Hfo")
surface_Hfo.setSpecificSurfaceArea(1e3, "m2/g").setMass(0.33, "g")

# Defined strong site of the complexation surface
surface_Hfo.addSite("Hfo_s", "_s").setAmount(1.0, "mol")

# Defined weak site of the complexation surface
site_Hfo_w = ComplexationSurfaceSite()
site_Hfo_w.setName("Hfo_w").setAmount(1.0, "mol")
surface_Hfo.addSite(site_Hfo_w)

# Add species to the surface and corresponding sites
surface_Hfo.addSurfaceSpecies(list)

# Add specified surface as parameters for the activity model for the complexation surface
params = ActivityModelSurfaceComplexationParams()
params.surface = surface_Hfo

# Define surface complexation phase and set an activity model
complexation_phase_Hfo = SurfaceComplexationPhase(list_str)
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationNoDDL(params))

print(surface_Hfo)


# Create chemical system
system = ChemicalSystem(db, solution, complexation_phase_Hfo)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define initial equilibrium state
state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
state.set("H2O"   , 1.00, "kg")
state.set("Sr+2"  , 1.00, "mmol")
state.set("Hfo_wOH"  , surface_Hfo.sites()["_w"].amount(), "mol")
state.set("Hfo_sOH"  , surface_Hfo.sites()["_s"].amount(), "mol")

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)
res = solver.solve(state)

print("*******************************************\n"
      "After equilibration:\n"
      "*******************************************")
print("succeed       = ", res.optima.succeeded)
print("solutionstate = \n", state)

aqprops = AqueousProps(state)
print("I  = %e mol/kgw" % float(aqprops.ionicStrength()))
print("pH = %f" % float(aqprops.pH()))
print("pE = %f" % float(aqprops.pE()))
print(aqprops)

surfprops = ComplexationSurfaceProps(surface_Hfo, state)
print(surfprops)
