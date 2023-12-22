from reaktoro import *

db = PhreeqcDatabase("phreeqc.dat")

solution = AqueousPhase(speciate("H2O BaCl2 Na2SO4"))
solution.set(ActivityModelPitzer())

minerals = MineralPhases("Barite")

system = ChemicalSystem(db, solution, minerals)

fluid1 = Material(system)
fluid1.add("H2O"  , 100.0, "g")
fluid1.add("BaCl2",  10.0, "g")
fluid1.add("Na2SO4",  1.0, "mg")

fluid2 = Material(system)
fluid2.add("H2O"   , 100.0, "g")
fluid2.add("Na2SO4",  10.0, "g")
fluid2.add("BaCl2",  1.0, "mg")

x = 0.0  # or when x = 1

mix = fluid1(1.0 - x, "kg") + fluid2(x, "kg")

options = EquilibriumOptions()
options.optima.output.active = True

state = mix.equilibrate(25, "°C", 1, "bar", options)

assert mix.result().succeeded(), f"Equilibration failed at x = {x}"
