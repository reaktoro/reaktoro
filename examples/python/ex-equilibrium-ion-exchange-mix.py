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
#   • Svetlana Kyas (23 November, 2021)
#
# and since revised by:
#   •
# -----------------------------------------------------------------------------

from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

# Define an aqueous phase
solution = AqueousPhase("H2O Na+ Cl- H+ OH- K+ Ca+2")
solution.setActivityModel(chain(
    ActivityModelHKF()
))

# Define an ion exchange phase
exchange = IonExchangePhase("CaX2 KX NaX")
exchange.setActivityModel(ActivityModelIonExchangeVanselow())

# Create chemical system
system = ChemicalSystem(db, solution, exchange)

T = 25.0 # temperature in celsius
P = 1.0  # pressure in bar

# Define initial equilibrium state
state = ChemicalState(system)
state.setTemperature(T, "celsius")
state.setPressure(P, "atm")

# To match the PHREEQC script:
# ---------------------------------------------------------------
# EXCHANGE 1
# CaX2 0.4065 # exchangeable Ca and K in mol
# KX 0.1871
# SOLUTION 1
# Na 1.2; Cl 1.2 charge # Na in solution, exchanges for K and Ca
# MIX 1; 1 1e-9 # ...Take 1e-9 of solution 1 = 1 μg water
# END
# ---------------------------------------------------------------

mix_scale = 1e-9
state.setSpeciesMass("H2O"   , 1.00 * mix_scale, "kg")
state.setSpeciesAmount("Na+" , 1.20 * mix_scale, "mmol")
state.setSpeciesAmount("Cl-" , 1.20 * mix_scale, "mmol")
# Define exchange species
state.setSpeciesAmount("KX"  , 0.1870, "mol")
state.setSpeciesAmount("CaX2", 0.4065, "mol")

# Define equilibrium solver and equilibrate given initial state with input conditions
solver = EquilibriumSolver(system)
solver.solve(state)

# Calculate aqueous and exchange properties
aqprops = AqueousProps(state)
exprops = IonExchangeProps(state)

print(aqprops)
print(exprops)

print("I  = %f mol/kgw" % aqprops.ionicStrength()[0])
print("pH = %f" % aqprops.pH()[0])
print("pE = %f\n" % aqprops.pE()[0])

print("n(CaX2) = %e mole" % exprops.speciesAmount("CaX2")[0])
print("n(KX)   = %e mole" % exprops.speciesAmount("KX")[0])
print("n(NaX)  = %e mole\n" % exprops.speciesAmount("NaX")[0])

print("m(Ca+2) = %e molal" % aqprops.speciesMolality("Ca+2")[0])
print("m(Cl-)  = %e molal" % aqprops.speciesMolality("Cl-")[0])
print("m(K+)   = %e molal" % aqprops.speciesMolality("K+")[0])
print("m(Na+)  = %e molal\n" % aqprops.speciesMolality("Na+")[0])
