# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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
#   • Svetlana Kyas (4 February 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------

import reaktoro as rkt

# Define Thermofun database
db = rkt.ThermoFunDatabase("cemdata18")

# Print out species in cemdata18-thermofun.json:
for species in db.species():
    print(species.name())

solution = rkt.AqueousPhase(rkt.speciate("H O K Na S Si Ca Mg Al C Cl"))

# Set up a and b parameters for the ionic species (KOH, b = 0.123, a = 3.67)
params = rkt.ActivityModelDebyeHuckelParams()
params.aiondefault = 3.67
params.biondefault = 0.123
params.bneutraldefault = 0.123

solution.setActivityModel(rkt.ActivityModelDebyeHuckel(params))

# Define minerals phases
minerals = rkt.MineralPhases("Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag "
                             "C2S C3A C3S C4AF Lim Gp Brc K2SO4 K2O Na2SO4 Na2O")

# Define solid phases
solidphase_C3AFS084H  = rkt.SolidPhase("C3FS0.84H4.32 C3AFS0.84H4.32") # AlFeSi-hydrogarnet_ss
solidphase_ettringite = rkt.SolidPhase("ettringite ettringite30") # Ettrignite_ss
solidphase_OH_SO4_AFm = rkt.SolidPhase("C4AH13 monosulphate12") # Monosulfate_ss
solidphase_CSHQ       = rkt.SolidPhase("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH") # CSH_ss

# Define chemical system by providing database, aqueous phase, minerals, and solid solutions
system = rkt.ChemicalSystem(db, solution, minerals,
                            solidphase_C3AFS084H,
                            solidphase_ettringite,
                            solidphase_OH_SO4_AFm,
                            solidphase_CSHQ)

# Specify conditions to be satisfied at chemical equilibrium
specs = rkt.EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

# Define equilibrium solver
solver = rkt.EquilibriumSolver(specs)

# Define initial equilibrium state
state = rkt.ChemicalState(system)
state.setSpeciesMass("H2O@", 40, "g") # water # water/binder = 0.4, 40g water per 100g of cement clinker
# clinker phases
state.setSpeciesMass("C2S" ,  9.70, "g") # belite
state.setSpeciesMass("C3A" ,  7.72, "g") # aluminate
state.setSpeciesMass("C3S" , 67.31, "g") # alite
state.setSpeciesMass("C4AF",  8.14, "g") # ferrite, (CaO)4(Al2O3)(Fe|3|2O3)
# additional
state.setSpeciesMass("Gp"    , 3.13, "g") # gypsum, CaSO4(H2O)2
state.setSpeciesMass("Cal"   , 0.10, "g") # calcite, CaCO3
state.setSpeciesMass("Lim"   , 0.93, "g") # lime, CaO
state.setSpeciesMass("Brc"   , 1.31, "g") # brucite, Mg(OH)2
state.setSpeciesMass("K2SO4" , 1.34, "g") # potasium-sulfate
state.setSpeciesMass("K2O"   , 0.05, "g") # potasium oxide
state.setSpeciesMass("Na2SO4", 0.21, "g") # sodium sulfate
state.setSpeciesMass("Na2O"  , 0.05, "g") # sodium oxide
state.setSpeciesMass("O2@"   , 0.15, "g") # oxygen to stabilize the system

# Define temperature and pressure
T = 20.0 # in Celsius
P = 1.0 # in bar

# Define conditions to be satisfied at chemical equilibrium
conditions = rkt.EquilibriumConditions(specs)
conditions.temperature(T, "celsius")
conditions.pressure(P, "bar")

# Equilibrate the initial state with given conditions and component amounts
res = solver.solve(state, conditions)
print("res (cemdata18) = ", res.optima.succeeded)

# Output the chemical state to a file
state.output("state-cemdata18.txt")

# Output chemical properties to a file
props = rkt.ChemicalProps (state)
props.output("props.txt")

# Output aqueous properties to a file
aprops = rkt.AqueousProps (state)
aprops.output("aprops.txt")

