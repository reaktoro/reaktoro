import numpy as np
import pytest
import reaktoro


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[1] / "data"


def test_convergence_with_liquid_phase(num_regression):
    temperature = 343.15
    pressure = 1323.9767e5

    db = reaktoro.Database(str(get_test_data_dir() / 'supcrt07_with_lumped_oil_full.xml'))

    editor = reaktoro.ChemicalEditor(db)

    aqueous_species = [
        "H2O(l)",
        "CO2(aq)",
        "HCO3-",
        "CO3--",
        "H2S(aq)",
        "HS-",
        "SO4--",
        "Ca++",
        "Na+",
        "Cl-",
        "Fe++",
        "Ba++",
        "Sr++",
    ]
    gaseous_species = ["CO2(g)", "H2S(g)", "CH4(g)"]
    oil_species = ["H2S(liq)", "CO2(liq)", "LumpOil(liq)"]
    minerals = ["Anhydrite", "Barite", "Calcite", "Celestite", "Siderite", "Pyrrhotite"]

    editor.addAqueousPhase(aqueous_species)
    editor.addGaseousPhase(gaseous_species)
    editor.addLiquidPhase(oil_species)
    for mineral in minerals:
        editor.addMineralPhase(mineral)

    system = reaktoro.ChemicalSystem(editor)

    problem = reaktoro.EquilibriumProblem(system)

    problem.setTemperature(temperature)
    problem.setPressure(pressure)

    state = reaktoro.ChemicalState(system)
    state.setTemperature(temperature)
    state.setPressure(pressure)

    # AQUEOUS
    state.setSpeciesAmount("H2O(l)", 1047530951.9776191)
    state.setSpeciesAmount("CO2(aq)", 20774.560623553392)
    state.setSpeciesAmount("HCO3-", 282912.10935910186)
    state.setSpeciesAmount("CO3--", 1483.0702563458142)
    state.setSpeciesAmount("H2S(aq)", 12166.740676127629)
    state.setSpeciesAmount("HS-", 41503.521359359736)
    state.setSpeciesAmount("SO4--", 0.00000000000000000)
    state.setSpeciesAmount("Ca++", 182556.98986975398)
    state.setSpeciesAmount("Na+", 8919454.7940101717)
    state.setSpeciesAmount("Cl-", 8955756.8340308685)
    state.setSpeciesAmount("Fe++", 0.00000000000000000)
    state.setSpeciesAmount("Ba++", 13319.520269138628)
    state.setSpeciesAmount("Sr++", 20875.710568363389)

    # MINERAL
    state.setSpeciesAmount("Calcite", 999131754.50533485)
    state.setSpeciesAmount("Anhydrite", 0.00000000000000000)
    state.setSpeciesAmount("Barite", 0.00000000000000000)
    state.setSpeciesAmount("Celestite", 0.00000000000000000)
    state.setSpeciesAmount("Siderite", 6732617546.7550077)
    state.setSpeciesAmount("Pyrrhotite", 0.00000000000000000)

    # GASEOUS
    state.setSpeciesAmount("CO2(g)", 574511.37111462699)
    state.setSpeciesAmount("H2S(g)", 13303.310614909715)
    state.setSpeciesAmount("CH4(g)", 700402204.48213911)

    # LIQUID
    state.setSpeciesAmount("CO2(liq)", 821887.77968393208)
    state.setSpeciesAmount("H2S(liq)", 240596.16100715694)
    state.setSpeciesAmount("LumpOil(liq)", 12676647027.415014)

    problem.addState(state)

    solver = reaktoro.EquilibriumSolver(system)

    state = reaktoro.ChemicalState(system)

    converged = solver.solve(state, problem).optimum.succeeded

    n = state.speciesAmounts()

    b = state.elementAmounts()

    output = {}
    output["Element amounts [mol]"] = np.asarray(b)
    output["Species amounts [mol]"] = n

    num_regression.check(output)

    assert converged
