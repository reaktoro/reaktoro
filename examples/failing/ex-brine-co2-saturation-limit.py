from reaktoro import *

db = PhreeqcDatabase("pitzer.dat")

system = ChemicalSystem(db,
    AqueousPhase(speciate("NaCl CO2 CaCO3")).set(ActivityModelPitzer()),
    GaseousPhase("CO2(g) H2O(g)").set(ActivityModelPengRobinsonPhreeqc()),  # using ActivityModelPengRobinsonPhreeqcOriginal works instead
    MineralPhases("Calcite")
)

nCO2 = 0.240246773970415467  # solubility of CO2 at the bubble limit the gas phase has 1e-6 moles

state = ChemicalState(system)
state.setTemperature(30.0, "°C")
state.setPressure(10, "bar")
state.set("H2O", 1.0, "kg")
state.set("Na+", 1.0, "mol")
state.set("Cl-", 1.0, "mol")
state.set("CO2(g)", nCO2, "mol")

solver = EquilibriumSolver(system)

result = solver.solve(state)

assert result.succeeded(), "Failed to solve equilibrium!"

nGas = state.speciesAmount("CO2(g)") + state.speciesAmount("H2O(g)")

print(f"Total number of moles of gas in the system is {nGas} mol but it should be 1e-6 mol")

# Total number of moles of gas in the system is 38.876001374389574 mol but
# should be 1e-6 mol.

# The ActivityModelPengRobinson is producing an artificial aqueous phase with
# predominant H2O, and this is becoming more stable (less Gibbs energy) than
# having the aqueous phase alone.
