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
import numpy as np
import math

# Define Thermofun database
db = ThermoFunDatabase("cemdata18")

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
ss_C3AFS084H.setName("ss_AlFeSi-hydrogarnet")
ss_ettringite.setName("ss_Ettrignite_ss")
ss_OH_SO4_AFm.setName("ss_Monosulfate_ss")
ss_CSHQ.setName("ss_C-S-H")

# Define chemical system by providing database, aqueous phase, minerals, and solid solutions
system = ChemicalSystem(db, solution, minerals, gaseous,
                        ss_C3AFS084H,
                        ss_ettringite,
                        ss_OH_SO4_AFm,
                        ss_CSHQ)
# # Print out species in cemdata18-thermofun.json:
# for phase in system.phases():
#     print(phase.name())
#
# input()

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

props = ChemicalProps(system)
aprops = AqueousProps(system)

opts = EquilibriumOptions()
opts.optima.output.active = False
opts.epsilon = 1e-13

# We define the materials for our equilibrium recipe
# Cement clinker composition from XRF as given in Lothenbach et al., (2008) recalculated for 100g
# Define water
water = Material(system)
water.add("H2O", 1000.0, "g")

mCaCO3 = 00.00
mCaO = 65.70

d_CaCO3 = 1.0
d_CO2   = 1.0

steps_num = 50

array_mCaCO3 = np.zeros(steps_num)

phase_ss_C3AFS084H  = math.nan * np.ones(steps_num)
phase_ss_ettringite = math.nan * np.ones(steps_num)
phase_ss_OH_SO4_AFm = math.nan * np.ones(steps_num)
phase_ss_CSHQ       = math.nan * np.ones(steps_num)
phase_Cal           = math.nan * np.ones(steps_num)
phase_hydrotalcite  = math.nan * np.ones(steps_num)
phase_Portlandite   = math.nan * np.ones(steps_num)
phase_hemicarbonate = math.nan * np.ones(steps_num)
phase_monocarbonate = math.nan * np.ones(steps_num)
phase_AmorSl        = math.nan * np.ones(steps_num)
phase_FeOOHmic      = math.nan * np.ones(steps_num)
phase_Gbs           = math.nan * np.ones(steps_num)
phase_Mag           = math.nan * np.ones(steps_num)

amount_C3FS084H432    = math.nan * np.ones(steps_num)
amount_C3AFS084H432   = math.nan * np.ones(steps_num)
amount_ettringite     = math.nan * np.ones(steps_num)
amount_ettringite30   = math.nan * np.ones(steps_num)
amount_C4AH13         = math.nan * np.ones(steps_num)
amount_monosulphate12 = math.nan * np.ones(steps_num)
amount_CSHQTobD       = math.nan * np.ones(steps_num)
amount_CSHQTobH       = math.nan * np.ones(steps_num)
amount_CSHQJenH       = math.nan * np.ones(steps_num)
amount_CSHQJenD       = math.nan * np.ones(steps_num)
amount_KSiOH          = math.nan * np.ones(steps_num)
amount_NaSiOH         = math.nan * np.ones(steps_num)

for i in range(0, steps_num):

    cement_clinker = Material(system)
    cement_clinker.add("SiO2" , 20.47, "g")
    cement_clinker.add("CaCO3", mCaCO3, "g")
    cement_clinker.add("CaO"  , mCaO, "g")
    cement_clinker.add("Al2O3",  4.90, "g")
    cement_clinker.add("Fe2O3",  3.20, "g")
    cement_clinker.add("K2O"  ,  0.79, "g")
    cement_clinker.add("Na2O" ,  0.42, "g")
    cement_clinker.add("MgO"  ,  1.80, "g")
    cement_clinker.add("SO3"  ,  2.29, "g")
    cement_clinker.add("CO2"  ,  0.26, "g")
    cement_clinker.add("O2"   ,  0.15, "g")

    # Define a cement mix of 0.5 water/binder
    cement_mix = Material(system)
    cement_mix = cement_clinker(100.0, "g") + water(50.0, "g")

    print(f"Calculations with m(CaCO3) = {mCaCO3} g and m(CaO) = {mCaO} g:")

    # Equilibrate cement mix
    state = cement_mix.equilibrate(20.0, "celsius", 1.0, "bar", opts)
    res = cement_mix.result()
    print("res (cemdata18, run with material) = ", res.optima.succeeded)
    #state.output("python-state-cemdata18_1.txt")

    if not res.optima.succeeded:
        # Define equilibrium solver
        solver = EquilibriumSolver(specs)

        # Equilibrate the resulting chemical state with equilibrium solver
        solver.setOptions(opts)
        res = solver.solve(state, conditions)
        print("res (cemdata18, run with solver) = ", res.optima.succeeded)
        #state.output("python-state-cemdata18_2.txt")

    # Update chemical and aqueous properties to a file
    props.update(state)
    aprops.update(state)

    if res.optima.succeeded:
        array_mCaCO3[i] = mCaCO3

        phase_ss_C3AFS084H[i] = props.phaseProps("ss_AlFeSi-hydrogarnet").amount()[0]
        phase_ss_ettringite[i] = props.phaseProps("ss_Ettrignite_ss").amount()[0]
        phase_ss_OH_SO4_AFm[i] = props.phaseProps("ss_Monosulfate_ss").amount()[0]
        phase_ss_CSHQ[i] = props.phaseProps("ss_C-S-H").amount()[0]

        phase_Cal[i] = props.phaseProps("Cal").amount()[0]
        phase_hydrotalcite[i] = props.phaseProps("hydrotalcite").amount()[0]
        phase_Portlandite[i] = props.phaseProps("Portlandite").amount()[0]
        phase_hemicarbonate[i] = props.phaseProps("hemicarbonate").amount()[0]
        phase_monocarbonate[i] = props.phaseProps("monocarbonate").amount()[0]
        phase_AmorSl[i] = props.phaseProps("Amor-Sl").amount()[0]
        phase_FeOOHmic[i] = props.phaseProps("FeOOHmic").amount()[0]
        phase_Gbs[i] = props.phaseProps("Gbs").amount()[0]
        phase_Mag[i] = props.phaseProps("Mag").amount()[0]

        amount_C3FS084H432[i] = state.speciesAmount("C3FS0.84H4.32")[0]
        amount_C3AFS084H432[i] = state.speciesAmount("C3AFS0.84H4.32")[0]
        amount_ettringite[i] = state.speciesAmount("ettringite")[0]
        amount_ettringite30[i] = state.speciesAmount("ettringite30")[0]
        amount_C4AH13[i] = state.speciesAmount("C4AH13")[0]
        amount_monosulphate12[i] = state.speciesAmount("monosulphate12")[0]
        amount_CSHQTobD[i] = state.speciesAmount("CSHQ-TobD")[0]
        amount_CSHQTobH[i] = state.speciesAmount("CSHQ-TobH")[0]
        amount_CSHQJenH[i] = state.speciesAmount("CSHQ-JenH")[0]
        amount_CSHQJenD[i] = state.speciesAmount("CSHQ-JenD")[0]
        amount_KSiOH[i] = state.speciesAmount("KSiOH")[0]
        amount_NaSiOH[i] = state.speciesAmount("NaSiOH")[0]

    #cement_clinker.add("CaCO3",   d_CaCO3, "g")
    #cement_clinker.add("CaO"  ,  -d_CO2  , "g")

    mCaCO3 += d_CaCO3
    mCaO -= d_CO2

print("phase_ss_C3AFS084H : ", phase_ss_C3AFS084H)
print("phase_ss_ettringite : ", phase_ss_ettringite)
print("phase_ss_OH_SO4_AFm : ", phase_ss_OH_SO4_AFm)
print("phase_ss_CSHQ : ", phase_ss_CSHQ)

print("phase_Cal : ", phase_Cal)
print("phase_hydrotalcite : ", phase_hydrotalcite)
print("phase_hemicarbonate : ", phase_hemicarbonate)
print("phase_monocarbonate : ", phase_monocarbonate)


import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9', 'C0', 'darkblue', 'darkgreen']

plt.figure()
plt.xlabel("Mass of CaCO3 [g]")
plt.ylabel("Phases amounts [mol]")

plt.plot(array_mCaCO3, phase_ss_C3AFS084H, label="ss_AlFeSi-hydrogarnet", color=colors[0])

plt.legend(loc="best")
plt.grid()
plt.savefig(f'phases-AlFeSi-hydrogarnet-vs-caco2.png', bbox_inches='tight')
plt.close()

# plt.figure()
# plt.xlabel("Mass of CaCO3 [g]")
# plt.ylabel("Phases amounts [mol]")
#
# plt.plot(array_mCaCO3, phase_ss_C3AFS084H, label="ss_AlFeSi-hydrogarnet", color=colors[0])
# plt.plot(array_mCaCO3, phase_ss_CSHQ, label="ss_C-S-H", color=colors[3])
#
# plt.legend(loc="best")
# plt.grid()
# plt.savefig(f'phases-AlFeSi-hydrogarnet-C-S-H-vs-caco2.png', bbox_inches='tight')
# plt.close()


# plt.figure()
# plt.xlabel("Mass of CaCO3 [g]")
# plt.ylabel("Phases amounts [mol]")
#
# plt.plot(array_mCaCO3, phase_ss_ettringite, label="ss_Ettrignite_ss", color=colors[1])
#
# plt.legend(loc="best")
# plt.grid()
# plt.savefig(f'phases-Ettrignite-vs-caco2.png', bbox_inches='tight')
# plt.close()
#
#
plt.figure()
plt.xlabel("Mass of CaCO3 [g]")
plt.ylabel("Phases amounts [mol]")

plt.plot(array_mCaCO3, phase_ss_CSHQ, label="ss_C-S-H", color=colors[3])

plt.legend(loc="best")
plt.grid()
plt.savefig(f'phases-CSH-vs-caco2.png', bbox_inches='tight')
plt.close()

plt.figure()
plt.xlabel("Mass of CaCO3 [g]")
plt.ylabel("Phases amounts [mol]")

plt.plot(array_mCaCO3, phase_Cal, label="Calcite", color=colors[2])

plt.legend(loc="best")
plt.grid()
#plt.show()
plt.savefig(f'phases-calcite-vs-caco2.png', bbox_inches='tight')
plt.close()


# plt.figure()
# plt.xlabel("Mass of CaCO3 [g]")
# plt.ylabel("Phases amounts [mol]")
#
# plt.plot(array_mCaCO3, phase_hydrotalcite, label="Hydrotalcite", color=colors[4])
#
# plt.legend(loc="best")
# plt.grid()
# #plt.show()
# plt.savefig(f'phases-hydrotalcite-vs-caco2.png', bbox_inches='tight')
# plt.close()
#
#
# plt.figure()
# plt.xlabel("Mass of CaCO3 [g]")
# plt.ylabel("Phases amounts [mol]")
#
# plt.plot(array_mCaCO3, phase_monocarbonate, label="Monocarbonate", color=colors[5])
#
# plt.legend(loc="best")
# plt.grid()
# #plt.show()
# plt.savefig(f'phases-monocarbonate-vs-caco2.png', bbox_inches='tight')
# plt.close()

plt.figure()
plt.xlabel("Mass of CaCO3 [g]")
plt.ylabel("Phases amounts [mol]")

plt.plot(array_mCaCO3, phase_monocarbonate, label="Monocarbonate", color=colors[5])
plt.plot(array_mCaCO3, phase_hydrotalcite, label="Hydrotalcite", color=colors[4])
plt.plot(array_mCaCO3, phase_ss_ettringite, label="ss_Ettrignite_ss", color=colors[1])

plt.legend(loc="best")
plt.grid()
#plt.show()
plt.savefig(f'phases-monocarbonate-hydrotalcite-Ettrignite-vs-caco2.png', bbox_inches='tight')
plt.close()
