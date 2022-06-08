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
solution = AqueousPhase(speciate("H O Cl Ca Sr"))
solution.setActivityModel(ActivityModelHKF())

#    Hfo_sOH           1.126e-05       0.450   1.126e-05      -4.949
#    Hfo_sOHCa+2       1.034e-05       0.413   1.034e-05      -4.986
#    Hfo_sOH2+         1.714e-06       0.069   1.714e-06      -5.766
#    Hfo_sO-           1.693e-06       0.068   1.693e-06      -5.771
#    Hfo_sOHSr+2       1.134e-11       0.000   1.134e-11     -10.945
#
#    Hfo_wOH           7.666e-04       0.767   7.666e-04      -3.115
#    Hfo_wOH2+         1.167e-04       0.117   1.167e-04      -3.933
#    Hfo_wO-           1.153e-04       0.115   1.153e-04      -3.938
#    Hfo_wOCa+         1.364e-06       0.001   1.364e-06      -5.865
#    Hfo_wOSr+         2.542e-13       0.000   2.542e-13     -12.595
#    Hfo_wOSrOH        3.108e-16       0.000   3.108e-16     -15.508

# Define ion exchange species list
list_str = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH"
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

print(surface_Hfo)

# Add specified surface as parameters for the activity model for the complexation surface
params = ActivityModelSurfaceComplexationParams()
params.surface = surface_Hfo

# Define surface complexation phase and set an activity model
complexation_phase_Hfo = SurfaceComplexationPhase(list_str)
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params))

# Define the DDL phase
ddl_Hfo = DoubleLayerPhase(speciate("H O Cl Ca Sr"))
ddl_Hfo.named("DoubleLayerPhase")
params_dll = ActivityModelDDLParams()
#dl_Hfo.setActivityModel(chain(ActivityModelHKF(), ActivityModelDDL(params_dll)))
ddl_Hfo.setActivityModel(ActivityModelDDL(params_dll))

# Create chemical system
system = ChemicalSystem(db, solution, complexation_phase_Hfo, ddl_Hfo)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define initial equilibrium state
state = ChemicalState(system)
state.temperature(T, "celsius")
state.pressure(P, "bar")
state.set("H2O" , 1.00, "kg")
state.set("Cl-" , 2e+0, "mmol")
state.set("Ca+2", 1e+0, "mmol")
state.set("Sr+2", 1e-6, "mmol")
state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

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
print("aqprops = \n", aqprops)

surfprops = ComplexationSurfaceProps(surface_Hfo, state)
print(surfprops)

dlprops = DoubleLayerProps(state)
print("dlprops = \n", dlprops)