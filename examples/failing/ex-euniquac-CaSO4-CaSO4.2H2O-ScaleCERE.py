import reaktoro as rkt

db = rkt.Database.fromFile("/home/allan/codes/reaktoro/ScaleCERE.yaml")

params = rkt.Params.local("/home/allan/codes/reaktoro/ScaleCERE.yaml")

solution = rkt.AqueousPhase()
solution.set(rkt.ActivityModelExtendedUNIQUAC(params))

minerals = rkt.MineralPhases("Anhydrite Gypsum")

system = rkt.ChemicalSystem(db, solution, minerals)

state = rkt.ChemicalState(system)
state.temperature(5.0, "°C")
state.pressure(1000.0, "bar")
state.set("H2O", 1.0, "kg")
state.set("Anhydrite", 10.0, "mol")

options = rkt.EquilibriumOptions()
options.optima.output.active = True

result = rkt.equilibrate(state, options)

print(f"The result should be True but currently is: **{result.succeeded()}**")
