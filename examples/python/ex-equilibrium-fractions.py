from reaktoro import *
import numpy as np

db = SupcrtDatabase('supcrtbl')
p = 1
t = 372.7559288971232 + 15
aqueous = AqueousPhase("H2O(aq)")
gaseous = GaseousPhase("H2O(g)")

system = ChemicalSystem(db, aqueous, gaseous)

opt = EquilibriumOptions()
opt.warmstart = False # default value is true
solver = EquilibriumSolver(system)
solver.setOptions(opt)

def initilized_state():
    state = ChemicalState(system)
    state.set('H2O(aq)', 1, 'kg')
    return state

def calc_rkt_frac(state, t, p):

    state.pressure(p, 'bar')
    state.temperature(t, 'K')

    try:
        res = solver.solve(state)
        print("Equilibration Succeeded? ", res.optima.succeeded)

        props = ChemicalProps(state)
        phases = state.system().phases().size()
        masses = np.array([props.phaseProps(i).mass()[0] for i in range(phases)])
        volumes = np.array([props.phaseProps(i).volume()[0] for i in range(phases)])

        densities = masses / (volumes + 1e-6)
        total_mass = sum(masses)

        fractions = masses / total_mass

        for i in range(len(fractions)):
            if densities[i] < 50:
                fractions[i] = 0.0
        fractions = fractions / (sum(fractions) + 1e-6)

        fraction = 1 - sum(fractions)

        return fraction

    except:
        print("RKT failed to converge with p={:.2f} and t={:.1f}".format(p, t))
        return -1


state = initilized_state()
print("Without re-building the state")
print(calc_rkt_frac(state, t, p))
print(calc_rkt_frac(state, t - 30, p))
print(calc_rkt_frac(state, t, p))

# recreate the state each time
print("\nRe-building the state each time")
state = initilized_state()
print(calc_rkt_frac(state, t, p))

state = initilized_state()
print(calc_rkt_frac(state, t - 30, p))

state = initilized_state()
print(calc_rkt_frac(state, t, p))
