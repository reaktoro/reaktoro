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
#   • G.D. Miron (1 April 2022)
# -----------------------------------------------------------------------------

import sys
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

# Define AlFeSi-hydrogarnet solid phase
ss_C3AFS084H  = SolidPhase("C3FS0.84H4.32 C3AFS0.84H4.32")
ss_C3AFS084H.setName("ss_C3AFS084H")
# Define Ettrignite solid phase
ss_ettringite = SolidPhase("ettringite ettringite30")
ss_ettringite.setName("ss_Ettrignite")
# Define Monosulfate solid phase
ss_OH_SO4_AFm = SolidPhase("C4AH13 monosulphate12")
ss_OH_SO4_AFm.setName("ss_Monosulfate")
# Define CSH solid phase
ss_CSHQ = SolidPhase("CSHQ-TobD CSHQ-TobH CSHQ-JenH CSHQ-JenD KSiOH NaSiOH")
ss_CSHQ.setName("ss_CSHQ")

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

# Define conditions to be satisfied at chemical equilibrium
conditions = EquilibriumConditions(specs)
conditions.temperature(20.0, "celsius")
conditions.pressure(1.0, "bar")

props = ChemicalProps(system)
aprops = AqueousProps(system)

opts = EquilibriumOptions()
opts.optima.output.active = False
opts.epsilon = 1e-13

# We define the materials for our equilibrium recipe
# Cement clinker composition from XRF as given in Lothenbach et al., (2008) recalculated for 100g
cement_clinker = Material(system)
cement_clinker.add("SiO2" , 20.47, "g")
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

# Define calcite
calcite = Material(system)
calcite.add("CaCO3", 1, "g")

import numpy as np
# Create list of species and phases names, list of Species objects, and auxiliary amounts array
phases_list_str = "ss_C3AFS084H ss_Ettrignite ss_Monosulfate ss_CSHQ " \
                  "Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag".split()
volume = np.zeros(len(phases_list_str))

# Define dataframe to collect amount of the selected species
import pandas as pd
columns = ["CaCO3"] \
          + ["volume_perc_" + name for name in phases_list_str]
df = pd.DataFrame(columns=columns)

# Volume in cm3
total_volume = float(props.volume()) *1e5

# Number of steps
steps_num = 19

import numpy as np
# Create list of species and phases names, list of Species objects, and auxiliary amounts array
phases_list_str = "ss_C3AFS084H ss_Ettrignite ss_Monosulfate ss_CSHQ " \
                  "Cal hydrotalcite Portlandite hemicarbonate monocarbonate Amor-Sl FeOOHmic Gbs Mag".split()
volume = np.zeros(len(phases_list_str))

# Define dataframe to collect amount of the selected species
import pandas as pd
columns = ["CaCO3"] \
          + ["volume_perc_" + name for name in phases_list_str]
df = pd.DataFrame(columns=columns)

# Number of steps
steps_num = 19

# We simulate the addition of calcite at the expense of clinker in the cement mix
for i in range(1, steps_num):

    # Define a cement mix of 0.5 water/binder at each step calcite is added at the expense of clinker
    cement_mix = Material(system)
    cement_mix = cement_clinker(100.0-i, "g") + calcite(i, "g") + water(50.0, "g")

    # Equilibrate cement mix
    state = cement_mix.equilibrate(20.0, "celsius", 1.0, "bar", opts)
    res = cement_mix.result()

    if not res.optima.succeeded:
        # Define equilibrium solver
        solver = EquilibriumSolver(specs)

        # Equilibrate the resulting chemical state with equilibrium solver
        solver.setOptions(opts)
        res = solver.solve(state, conditions)

        if not res.optima.succeeded: continue

    # Update chemical and aqueous properties to a file
    props.update(state)
    aprops.update(state)

    for j in range(0, len(phases_list_str)):
        # Collecting volume in cm3
        volume[j] = float(props.phaseProps(phases_list_str[j]).volume())
    volume_perc = volume / float(props.volume()) * 100

    # Update dataframe with obtained values
    df.loc[len(df)] = np.concatenate([[i], volume_perc])

# Plot selected dataframe columns
volume_names = ["Cal", "hydrotalcite", "Portlandite", "ss_CSHQ", "ss_C3AFS084H", "ss_Ettrignite", "monocarbonate"]
ax = df.plot.area(x='CaCO3', y=["volume_perc_" + name for name in volume_names], label=volume_names, colormap="Set2")
ax.set_xlim(2.0, 10.0)
ax.figure.savefig('phase-volume-vs-caco3.png')
