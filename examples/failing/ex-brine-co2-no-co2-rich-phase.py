from reaktoro import *

db = PhreeqcDatabase("pitzer.dat")

system = ChemicalSystem(db,
    AqueousPhase(speciate("NaCl CO2 CaCO3")).set(ActivityModelPitzer()),
    GaseousPhase("CO2(g) H2O(g)").set(ActivityModelPengRobinsonPhreeqc()),  # using ActivityModelPengRobinsonPhreeqcOriginal works instead
    MineralPhases("Calcite")
)

nCO2 = 1.0

state = ChemicalState(system)
state.setTemperature(30.0, "°C")
state.setPressure(10, "bar")
state.set("H2O", 1.0, "kg")
state.set("Na+", 1.0, "mol")
state.set("Cl-", 1.0, "mol")
state.set("CO2", nCO2, "mol")

# Uncomment these lines and the expected result below is obtained
# state.set("CO2(g)", 1e-14, "mol")
# state.set("H2O(g)", 1e-15, "mol")

solver = EquilibriumSolver(system)

result = solver.solve(state)

assert result.succeeded(), "Failed to solve equilibrium!"

nGas = state.speciesAmount("CO2(g)") + state.speciesAmount("H2O(g)")

print(f"Total number of moles of gas in the system is {nGas} mol but it should be 1e-6 mol.")

# Total number of moles of gas in the system is 2e-16 mol but should be
# 0.7631740905950601 mol. The CO2-rich phase is not forming.
# ActivityModelPengRobinson is functioning as if the gaseous phase was liquid,
# because the initial condition of the gases are 1e-16 mol each, and thus 0.5
# mole fraction.
