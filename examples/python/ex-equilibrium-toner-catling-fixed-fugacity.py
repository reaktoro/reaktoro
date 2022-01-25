# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: notebooks//ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
# ---

from reaktoro import *
import numpy as np
import os
import math

results_folder = 'results-phrqc2-fugacity-fixed-different-ppCO2'
os.system('mkdir -p ' + results_folder)

db = PhreeqcDatabase.fromFile('../examples/resources/phreeqc-toner-catling.dat') # if running from tutorials folder

solution = AqueousPhase(speciate("H O C Na Cl"))
solution.setActivityModel(chain(
    ActivityModelHKF(),
    ActivityModelDrummond("CO2")
))

minerals = MineralPhases("Natron Nahcolite Trona Na2CO3:H2O Na2CO3:7H2O")

system = ChemicalSystem(db, solution, minerals)

props = ChemicalProps(system)
aprops = AqueousProps(system)

specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.fugacity("CO2")

conditions = EquilibriumConditions(specs)

T = 50.0 # temperature in celsius
P = 1.0  # pressure in bar

solver = EquilibriumSolver(specs)

opts = EquilibriumOptions()
opts.epsilon = 1e-13
solver.setOptions(opts)

state0 = ChemicalState(system)
state0.set("H2O"        ,   1.00, "kg")
state0.set("Nahcolite"  ,  10.00, "mol")    # NaHCO3
state0.set("Natron"     ,   0.00, "mol")    # Na2CO3:10H2O
state0.set("Trona"      ,   0.00, "mol")    # Na3H(CO3)2:2H2O
state0.set("Na2CO3:H2O" ,   0.00, "mol")
state0.set("Na2CO3:7H2O",   0.00, "mol")
state0.set("CO2"        , 100.00, "mol")

# Explanation of the issues:
num_ppco2s = 17; initialize = False # everything is converging for all the temperatures (if state is NOT initialized in the loop)
#num_ppco2s = 17; initialize = True # ppCO2 = -2.375 does not converge (if state is initialized in the loop)
#num_ppco2s = 71; initialize = False # starting from ppCO2 = -1.9, does not converge for T = 50C
#num_ppco2s = 71; initialize = True # BUT! only for ppCO2 = -2.3, does not converge for T = 50C (if state is initialized in the loop)

co2ppressures = np.flip(np.linspace(-5.0, 2.0, num=num_ppco2s))

data_size = 3
data = np.zeros((num_ppco2s, data_size+1))

def equilibrate(ppCO2, state):

    conditions.temperature(T, "celsius")
    conditions.pressure(P, "atm")
    conditions.fugacity("CO2", 10**(ppCO2), "atm")

    if initialize:
        state = ChemicalState(system)
        state.set("H2O"        ,   1.00, "kg")
        state.set("Nahcolite"  ,  10.00, "mol")    # NaHCO3
        state.set("Natron"     ,   0.00, "mol")    # Na2CO3:10H2O
        state.set("Trona"      ,   0.00, "mol")    # Na3H(CO3)2:2H2O
        state.set("Na2CO3:H2O" ,   0.00, "mol")
        state.set("Na2CO3:7H2O",   0.00, "mol")
        state.set("CO2"        , 100.00, "mol")

    res = solver.solve(state, conditions)

    if not res.optima.succeeded:
        print(f"The optimization solver hasn't converged for T = {T} C and ppCO2 = {ppCO2}")
        return math.nan, math.nan, math.nan

    props.update(state)
    aprops.update(state)

    ph = aprops.pH()[0]
    mCO3 = state.speciesAmount("CO3-2")[0]
    mHCO3 = state.speciesAmount("HCO3-")[0]

    return ph, mCO3, mHCO3

print(f' ppCO2     pH      CO3-2      HCO3-')
for i in range(0, num_ppco2s):
    state = state0
    result = equilibrate(co2ppressures[i], state)
    data[i, 0] = co2ppressures[i]
    data[i, 1] = result[0]
    data[i, 2] = result[1]
    data[i, 3] = result[2]
    print(f'{co2ppressures[i]:8.4f} {result[0]:6.2f} {result[1]:6.4e} {result[2]:6.4e}')

np.savetxt(results_folder + '/m-data.txt', data)
np.savetxt(results_folder + '/m-pH.txt', data[:, 1])

import matplotlib.pyplot as plt
colors = ['C1', 'C2', 'C3', 'C4', 'C5', 'C7', 'C8', 'C9']

plt.figure()
plt.plot(data[:, 0], data[:, 1], label=f'pH', color=colors[1])
plt.legend(loc="best")
plt.xlabel('ppCO2')
plt.ylabel('pH [-]')
plt.grid()
plt.savefig(results_folder + '/' + 'pH-vs-ppCO2.png', bbox_inches='tight')
plt.close()