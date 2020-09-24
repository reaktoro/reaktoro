# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2020 Allan Leal
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


from reaktoro import *

import numpy as np
import matplotlib.pyplot as plt


def solubility_co2(system, T, P, mNaCl):

    n0CO2g = 10.0

    state = ChemicalState(system)

    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesMass("H2O", 1.0, "kg")
    state.setSpeciesAmount("CO2(g)", n0CO2g, "mol")
    state.setSpeciesAmount("Na+", mNaCl, "mol")
    state.setSpeciesAmount("Cl-", mNaCl, "mol")

    solver = EquilibriumSolver(system)

    solver.solve(state)

    nCO2g = state.speciesAmount("CO2(g)")

    return n0CO2g - nCO2g



db = PhreeqcDatabase("phreeqc.dat")

aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
aqueousphase.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2"),
))

gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.setActivityModel(ActivityModelPengRobinson())

phases = Phases(db)
phases.add(aqueousphase)
phases.add(gaseousphase)

system = ChemicalSystem(phases)

T = np.arange(25.0, 90.0, 5.0)
P = 1.0

mCO2_1 = [solubility_co2(system, x, P, mNaCl=1.0)[0] for x in T]  # [0] is needed to get the value of autodiff.real
mCO2_2 = [solubility_co2(system, x, P, mNaCl=2.0)[0] for x in T]  # [0] is needed to get the value of autodiff.real
mCO2_3 = [solubility_co2(system, x, P, mNaCl=4.0)[0] for x in T]  # [0] is needed to get the value of autodiff.real

fig, ax = plt.subplots()

ax.plot(T, mCO2_1, label=f"1 NaCl molal")
ax.plot(T, mCO2_2, label=f"2 NaCl molal")
ax.plot(T, mCO2_3, label=f"4 NaCl molal")

ax.legend(loc="upper right")

ax.set(xlabel='Temperature [degC]', ylabel='Solubility [mol/kgw]',
       title='Solubility of CO2 in NaCl brine')
ax.grid()

fig.savefig("co2-solubility-nacl-h2o.png")
plt.show()
