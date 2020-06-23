from collections import namedtuple

SpeciesAmount = namedtuple('SpeciesAmount', ('name', 'index', 'molar_amount'))
ReaktoroCase = namedtuple('ReaktoroCase', ('pressure_in_Pa', 'temperature_in_K', 'species_amounts'))


def get_reaktoro_case():
    return ReaktoroCase(
        pressure_in_Pa=13660337.317389656,
        temperature_in_K=64.0 + 273.15,
        species_amounts=(
            SpeciesAmount("H2O(l)", 32, 22360.14714403098),
            SpeciesAmount("CO2(aq)", 6, 25.350676739651337),
            SpeciesAmount("HCO3-", 37, 40.700744030286614),
            SpeciesAmount("CO3--", 7, 0.0265024560320966),
            SpeciesAmount("H2S(aq)", 34, 0),
            SpeciesAmount("HS-", 44, 0),
            SpeciesAmount("SO4--", 70, 0.062485989848723215),
            SpeciesAmount("Ca++", 9, 0.27957460644472432),
            SpeciesAmount("Na+", 50, 135.22073966783472),
            SpeciesAmount("Cl-", 15, 94.900737455906253),
            SpeciesAmount("Fe++", 20, 0.0041754804501114755),
            SpeciesAmount("Ba++", 1, 0.058280033772838649),
            SpeciesAmount("Sr++", 72, 0.20095330786840615),
            SpeciesAmount("Calcite", 81, 999131754.50533485),
            SpeciesAmount("Anhydrite", 79, 0),
            SpeciesAmount("Barite", 80, 0),
            SpeciesAmount("Celestite", 82, 0),
            SpeciesAmount("Siderite", 83, 0.01795172676219817),
            SpeciesAmount("Pyrrhotite", 84, 0),
        ),
    )


def test_different_results(state_regression):
    from reaktoro import ChemicalEditor, ChemicalState, ChemicalSystem, Database, EquilibriumProblem, EquilibriumOptions, EquilibriumSolver, Partition

    database = Database('supcrt07.xml')
    editor = ChemicalEditor(database)

    aqueous_elements = ["C", "Ca", "Cl", "Fe", "H", "Na", "O", "S", "Ba", "Sr"]
    aqueous_phase = editor.addAqueousPhaseWithElements(aqueous_elements)
    assert aqueous_phase.name() == 'Aqueous'

    mineral_species = ["Anhydrite", "Barite", "Calcite", "Celestite", "Siderite", "Pyrrhotite"]
    for mineral in mineral_species:
        editor.addMineralPhase(mineral)

    gaseous_species = ["CO2(g)", "H2S(g)", "CH4(g)"]
    editor.addGaseousPhase(gaseous_species)

    chemical_system = ChemicalSystem(editor)

    element_index = {e.name(): index for index, e in enumerate(chemical_system.elements())}
    species_index = {s.name(): index for index, s in enumerate(chemical_system.species())}
    phase_index = {p.name(): index for index, p in enumerate(chemical_system.phases())}

    reaktoro_case = get_reaktoro_case()

    partition = Partition(chemical_system)
    partition.setInertPhases([phase_index['Gaseous']])

    equilibrium_problem = EquilibriumProblem(partition)
    equilibrium_problem.setTemperature(reaktoro_case.temperature_in_K)
    equilibrium_problem.setPressure(reaktoro_case.pressure_in_Pa)

    chemical_state = ChemicalState(chemical_system)
    for name, index, molar_amount in reaktoro_case.species_amounts:
        assert index == species_index[name]
        chemical_state.setSpeciesAmount(index, molar_amount)

    equilibrium_problem.addState(chemical_state)

    options = EquilibriumOptions()
    options.optimum.output.active = True

    solver = EquilibriumSolver(partition)
    solver.setOptions(options)

    result = solver.solve(chemical_state, equilibrium_problem)

    assert result.optimum.succeeded

    state_regression.check(chemical_state, default_tol=dict(atol=1e-5, rtol=1e-14))
