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
complexation_phase_Hfo.setActivityModel(ActivityModelSurfaceComplexationWithDDL(params))

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

# Define initial equilibrium state
state = ChemicalState(system)
state.set("H2O" , 1.00, "kg")
state.set("Cl-" , 2e+0, "mmol")
state.set("Ca+2", 1e+0, "mmol")
state.set("Cd+2", 1e-6, "mmol")
state.set("Hfo_wOH", surface_Hfo.sites()["_w"].amount(), "mol")
state.set("Hfo_sOH", surface_Hfo.sites()["_s"].amount(), "mol")

# Equilibrate given initial state with input conditions
res = solver.solve(state, conditions)

# Update properties
props.update(state)

mCd_aq   = float(props.elementAmountInPhase("Cd", "AqueousPhase"))
mCd_surf = float(props.elementAmountInPhase("Cd", "SurfaceComplexationPhase"))
mCd_tot  = float(props.elementAmount("Cd"))

print(f"Total(Cd)     = {mCd_tot:6.4e} \n"
      f"Sorbed(Cd)    = {mCd_surf:6.4e}\n"
      f"Dissolved(Cd) = {mCd_aq:6.4e} \n "
      f"Sorbed(Cd)/Dissolved(Cd) = {mCd_surf/mCd_aq:6.4e} \n "
      f"Sorbed(Cd)/Total(Cd)     = {mCd_surf/mCd_tot:6.4e} ")