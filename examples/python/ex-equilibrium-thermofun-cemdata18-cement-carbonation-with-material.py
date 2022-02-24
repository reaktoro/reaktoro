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
#   • Svetlana Kyas (24 February 2022)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------

from reaktoro import *

# Define Thermofun database
db = ThermoFunDatabase("cemdata18")

# Print out species in cemdata18-thermofun.json:
for species in db.species():
    print(species.name())

solution = AqueousPhase(speciate("H O K Na S Si Ca Mg Al C Cl"))

# Set up a and b parameters for the ionic species (KOH, b = 0.123, a = 3.67)
params = ActivityModelDebyeHuckelParams()
params.aiondefault = 3.67
params.biondefault = 0.123
params.bneutraldefault = 0.123

solution.setActivityModel(ActivityModelDebyeHuckel(params))

# Define gas phase
gaseous = GaseousPhase(speciate("H O C"))

# Define minerals phases
minerals = MineralPhases("Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag ")

# Define solid phases
ss_C3AFS084H  = SolidPhase("C3FS0.84H4.32 C3AFS0.84H4.32") # AlFeSi-hydrogarnet_ss
ss_ettringite = SolidPhase("ettringite ettringite30") # Ettrignite_ss
ss_OH_SO4_AFm = SolidPhase("C4AH13 monosulphate12") # Monosulfate_ss
ss_CSHQ       = SolidPhase("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH") # CSH_ss

# Define chemical system by providing database, aqueous phase, minerals, and solid solutions
system = ChemicalSystem(db, solution, minerals, gaseous,
                        ss_C3AFS084H,
                        ss_ettringite,
                        ss_OH_SO4_AFm,
                        ss_CSHQ)

# Specify conditions to be satisfied at chemical equilibrium
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()

# Define temperature and pressure
T = 20.0 # in Celsius
P = 1.0 # in bar

# Define conditions to be satisfied at chemical equilibrium
conditions = EquilibriumConditions(specs)
conditions.temperature(T, "celsius")
conditions.pressure(P, "bar")

# Define equilibrium solver
solver = EquilibriumSolver(specs)
opts = EquilibriumOptions()
opts.optima.output.active = False

# We define the materials for our equilibrium recipe
# Cement clinker composition from XRF as given in Lothenbach et al., (2008) recalculated for 100g
cement_clinker = Material(system)
cement_clinker.add("SiO2" , 20.47, "g")
cement_clinker.add("CaCO3", 00.00, "g")
cement_clinker.add("CaO"  , 65.70, "g")
cement_clinker.add("Al2O3",  4.90, "g")
cement_clinker.add("Fe2O3",  3.20, "g")
cement_clinker.add("K2O"  ,  0.79, "g")
cement_clinker.add("Na2O" ,  0.42, "g")
cement_clinker.add("MgO"  ,  1.80, "g")
cement_clinker.add("SO3"  ,  2.29, "g")
cement_clinker.add("CO2"  ,  0.26, "g")
cement_clinker.add("O2"   ,  0.15, "g")

# Define water
water = Material(system)
water.add("H2O", 1000.0, "g")

# Define a cement mix of 0.5 water/binder
cement_mix = Material(system)
cement_mix = cement_clinker(100.0, "g") + water(50.0, "g")

# Equilibrate cement mix
state = cement_mix.equilibrate(20.0, "celsius", 1.0, "bar", opts)
res = cement_mix.result()
print("res (cemdata18, run with material) = ", res.optima.succeeded)
state.output("python-state-cemdata18_1.txt")

# Equilibrate the resulting chemical state with equilibrium solver
solver.setOptions(opts)
res = solver.solve(state, conditions)
print("res (cemdata18, run with solver) = ", res.optima.succeeded)
state.output("python-state-cemdata18_2.txt")

# Output chemical properties to a file
props = ChemicalProps (state)
props.output("props.txt")

# Output aqueous properties to a file
aprops = AqueousProps (state)
aprops.output("aprops.txt")

