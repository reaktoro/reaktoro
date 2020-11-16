from reaktoro import *

db = Database("supcrt98.xml")

editor = ChemicalEditor(db)
aqueousphase = editor.addAqueousPhaseWithElements("H O Na Cl")
aqueousphase.setActivityModelSetschenow("NaCl(aq)", 1.0)
editor.addMineralPhase("Halite")

system = ChemicalSystem(editor)

problem = EquilibriumProblem(system)
problem.setTemperature(30.0, "celsius")
problem.setPressure(300.0, "bar")
problem.add("H2O", 1.0, "kg")
problem.add("NaCl", 14.0, "mol")

state = equilibrate(problem)

print(state)