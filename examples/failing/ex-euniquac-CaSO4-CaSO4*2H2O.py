import reaktoro as rkt

db = rkt.Database.fromFile("euniquac-GarciaThesis2005.yaml")

params = rkt.Params.embedded("ExtendedUNIQUAC/GarciaThesis2005.yaml")

solution = rkt.AqueousPhase()
solution.set(rkt.ActivityModelExtendedUNIQUAC(params))

minerals = rkt.MineralPhases("CaSO4 CaSO4*2H2O")

system = rkt.ChemicalSystem(db, solution, minerals)

state = rkt.ChemicalState(system)
state.temperature(20.0, "°C")
state.pressure(1000.0, "bar")
state.set("H2O", 1.0, "kg")
state.set("CaSO4", 10.0, "mol")

result = rkt.equilibrate(state)

print(f"The result should be True but currently is: **{result.succeeded()}**")

state = rkt.ChemicalState(system)
state.temperature(15.0, "°C")
state.pressure(10.0, "bar")
state.set("H2O", 1.0, "kg")
state.set("CaSO4", 10.0, "mol")

result = rkt.equilibrate(state)

print(f"The result should be True but currently is: **{result.succeeded()}**")

state = rkt.ChemicalState(system)
state.temperature(30.0, "°C")
state.pressure(10.0, "bar")
state.set("H2O", 1.0, "kg")
state.set("CaSO4", 10.0, "mol")

result = rkt.equilibrate(state)

print(f"This calculation is working however: **{result.succeeded()}**")
