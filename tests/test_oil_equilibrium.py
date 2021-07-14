import numpy as np
import pytest
import reaktoro


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[0] / "data"


def _add_hydrocarbons_to_database(db):
    hydrocarbons_db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

    element_names_in_db = set(e.name() for e in db.elements())
    for element in hydrocarbons_db.elements():
        if element.name() not in element_names_in_db:
            db.addElement(element)

    for species in hydrocarbons_db.gaseousSpecies():
        db.addGaseousSpecies(species)

    for species in hydrocarbons_db.liquidSpecies():
        db.addLiquidSpecies(species)


def test_ternary_c1_c4_c10_bips_setup():
    """
    This is a ternary example from Whitson monograph. Retrieved from Problem 18 in
    Appendix B. This test compares the results defining BIPs as 'zero' matrix with
    the result without BIPs setup. Such results should be the same.

    .. reference:
        Whitson, C. H., & Brulé, M. R. (2000). Phase behavior (Vol. 20).
        Richardson, TX: Henry L. Doherty Memorial Fund of AIME, Society of Petroleum Engineers.
    """
    temperature = 280  # degF
    pressure = 500  # psi

    composition = np.array([0.5, 0.42, 0.08])  # molar fractions

    db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

    editor = reaktoro.ChemicalEditor(db)
    editor_bips = reaktoro.ChemicalEditor(db)

    gaseous_species = ["C1(g)", "C4(g)", "C10(g)"]
    oil_species = ["C1(liq)", "C4(liq)", "C10(liq)"]

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
        phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
    )

    editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

    def calculate_bips(T):
        k = np.zeros((3, 3))
        bips = reaktoro.BinaryInteractionParams(k)
        return bips

    eos_params_bips = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
        phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
        binary_interaction_values=calculate_bips
    )

    editor_bips.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params_bips)
    editor_bips.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params_bips)

    system = reaktoro.ChemicalSystem(editor)
    system_bips = reaktoro.ChemicalSystem(editor_bips)

    problem = reaktoro.EquilibriumProblem(system)
    problem_bips = reaktoro.EquilibriumProblem(system_bips)

    problem.setTemperature(temperature, 'degF')
    problem.setPressure(pressure, 'psi')
    problem.setElementAmount('C1', composition[0])
    problem.setElementAmount('C4', composition[1])
    problem.setElementAmount('C10', composition[2])

    problem_bips.setTemperature(temperature, 'degF')
    problem_bips.setPressure(pressure, 'psi')
    problem_bips.setElementAmount('C1', composition[0])
    problem_bips.setElementAmount('C4', composition[1])
    problem_bips.setElementAmount('C10', composition[2])

    # Custom initial guess below. This is necessary for hydrocarbon mixtures.
    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)  # start will all having 0.001 moles
    state.setSpeciesAmount("C1(g)", composition[0])  # overwrite amount of C1(g) (same below)
    state.setSpeciesAmount("C4(g)", composition[1])  
    state.setSpeciesAmount("C10(g)", composition[2])

    state_bips = reaktoro.ChemicalState(system_bips)
    state_bips.setSpeciesAmounts(0.001)
    state_bips.setSpeciesAmount("C1(g)", composition[0])
    state_bips.setSpeciesAmount("C4(g)", composition[1])  
    state_bips.setSpeciesAmount("C10(g)", composition[2])

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact

    solver = reaktoro.EquilibriumSolver(system)
    solver.setOptions(options)
    result = solver.solve(state, problem)

    solver_bips = reaktoro.EquilibriumSolver(system_bips)
    solver_bips.setOptions(options)
    result_bips = solver_bips.solve(state_bips, problem_bips)
    
    assert result.optimum.succeeded
    assert result_bips.optimum.succeeded

    # Compare results with and without setting BIPs
    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')
    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions = np.array([gas_phase_molar_fraction, liquid_phase_molar_fraction])

    molar_base_bips = state_bips.phaseAmount('Gaseous') + state_bips.phaseAmount('Liquid')
    gas_phase_molar_fraction_bips = state_bips.phaseAmount('Gaseous') / molar_base_bips
    liquid_phase_molar_fraction_bips = state_bips.phaseAmount('Liquid') / molar_base_bips
    phase_fractions_bips = np.array([gas_phase_molar_fraction_bips, liquid_phase_molar_fraction_bips])
    assert phase_fractions == pytest.approx(phase_fractions_bips)


def test_error_bips_setup():
    """
    The goal of this test is to check a BIP input with incompatible dimensions.
    The BIPs compose a symmetric matrix of size (n_species, n_species). This test checks
    for the case where the k BIPs matrix is not symmetric.
    """
    db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

    editor_bips = reaktoro.ChemicalEditor(db)

    gaseous_species = ["C1(g)", "C4(g)", "C10(g)"]

    def calculate_bips(T):
        one = 1.0
        # Wrong BIP matrix
        k = [
            [one, one, one],
            [one, one, one],
            [-one, one, one],
        ]
        bips = reaktoro.BinaryInteractionParams(k)
        return bips

    eos_params_bips = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
        phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
        binary_interaction_values=calculate_bips
    )
    with pytest.raises(RuntimeError, match='k is not symmetric'):
        editor_bips.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params_bips)


@pytest.mark.parametrize(
    "P, T, F_expected, x_expected, y_expected",
    [
        # psi, degF, [mol / mol], [mol / mol], [mol / mol]
        (500, 280, [0.853401, 1 - 0.853401], [0.08588, 0.46349, 0.45064], [0.57114, 0.41253, 0.01633]),
        (1500, 280, [0.566844, 1 - 0.566844], [0.330082, 0.513307, 0.156611], [0.629843, 0.348699, 0.021457]),
    ]
)
def test_ternary_c1_c4_c10_mixture(P, T, F_expected, x_expected, y_expected):
    """
    This is a ternary example from Whitson monograph. Retrieved from Problem 18 in
    Appendix B.

    .. reference:
        Whitson, C. H., & Brulé, M. R. (2000). Phase behavior (Vol. 20).
        Richardson, TX: Henry L. Doherty Memorial Fund of AIME, Society of Petroleum Engineers.
    """
    temperature = T  # degF
    pressure = P  # psi

    composition = np.array([0.5, 0.42, 0.08])  # molar fractions

    db = reaktoro.Database('supcrt07.xml')
    _add_hydrocarbons_to_database(db)

    editor = reaktoro.ChemicalEditor(db)

    gaseous_species = ["C1(g)", "C4(g)", "C10(g)"]
    oil_species = ["C1(liq)", "C4(liq)", "C10(liq)"]

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
    )

    editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    problem = reaktoro.EquilibriumProblem(system)

    problem.setTemperature(temperature, 'degF')
    problem.setPressure(pressure, 'psi')
    problem.setElementAmount('C1', composition[0])
    problem.setElementAmount('C4', composition[1])
    problem.setElementAmount('C10', composition[2])

    # Custom initial guess below. This is necessary for hydrocarbon mixtures.
    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)  # start will all having 0.001 moles
    state.setSpeciesAmount("C1(g)", composition[0])  # overwrite amount of C1(g) (same below)
    state.setSpeciesAmount("C4(g)", composition[1])  
    state.setSpeciesAmount("C10(g)", composition[2])

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact

    solver = reaktoro.EquilibriumSolver(system)
    solver.setOptions(options)
    result = solver.solve(state, problem)
    
    assert result.optimum.succeeded

    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')
    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    liquid_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions = np.array([gas_phase_molar_fraction, liquid_molar_fraction])
    assert phase_fractions == pytest.approx(F_expected, rel=5e-3)

    c1_gas_fraction = state.speciesAmount('C1(g)') / phase_fractions[0]
    c4_gas_fraction = state.speciesAmount('C4(g)') / phase_fractions[0]
    c10_gas_fraction = state.speciesAmount('C10(g)') / phase_fractions[0]
    equilibrium_gas_composition = np.array([
        c1_gas_fraction,
        c4_gas_fraction,
        c10_gas_fraction
    ])
    assert equilibrium_gas_composition == pytest.approx(
        y_expected,
        rel=2e-2
    )

    c1_liq_fraction = state.speciesAmount('C1(liq)') / phase_fractions[1]
    c4_liq_fraction = state.speciesAmount('C4(liq)') / phase_fractions[1]
    c10_liq_fraction = state.speciesAmount('C10(liq)') / phase_fractions[1]
    equilibrium_liq_composition = np.array([
        c1_liq_fraction,
        c4_liq_fraction,
        c10_liq_fraction
    ])
    assert equilibrium_liq_composition == pytest.approx(
        x_expected,
        rel=2e-2
    )


def test_equilibrium_CH4_CO2_pedersen():
    """
    This case is gathered from Pedersen book [1]. It is a C1-CO2 mixture using
    SRK as Cubic EOS. In P = 20 bar and T = -42 degC, it results in a LV equilibrium.
    This test compare the results from the book with Reaktoro's outcomes. The main
    goal is to check if BIPs setup is working properly. The case is given in Table 6.3
    in the Chapter 6 of [1]. The BIP value is given in the Chapter 4, Table 4.2.

    .. reference:
        [1] Pedersen, K. S., Christensen, P. L., & Shaikh, J. A. (2014). Phase behavior 
            of petroleum reservoir fluids. CRC press.
    """
    
    db = reaktoro.Database("supcrt98.xml") 
     
    editor = reaktoro.ChemicalEditor(db)

    def calculate_bips(T):
        k = np.zeros((2, 2))
        k[0, 1] = k[1, 0] = 0.12
        bips = reaktoro.BinaryInteractionParams(k)
        return bips

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.SoaveRedlichKwong,
        phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
        binary_interaction_values=calculate_bips
    )
    
    editor.addGaseousPhase(["CH4(g)", "CO2(g)"]).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(["CH4(liq)", "CO2(liq)"]).setChemicalModelCubicEOS(eos_params)
       
    system = reaktoro.ChemicalSystem(editor)
    
    problem = reaktoro.EquilibriumProblem(system)
    
    temperature = -42.0  # degC
    pressure = 20.0  # bar
    problem.setTemperature(temperature, "degC")
    problem.setPressure(pressure, "bar")
    problem.add("CH4(g)", 0.4, "mol")
    problem.add("CO2(g)", 0.6, "mol")

    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)  # start will all having 0.001 moles
    state.setSpeciesAmount("CH4(g)", 0.4)  # overwrite amount of C1(g) (same below)
    state.setSpeciesAmount("CO2(g)", 0.6)
    
    solver = reaktoro.EquilibriumSolver(problem.system())
    
    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact
    solver.setOptions(options)
    
    result = solver.solve(state, problem)

    assert result.optimum.succeeded

    # Compare phases molar fractions
    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')
    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions = np.array([liquid_phase_molar_fraction, gas_phase_molar_fraction])
    F_expected = np.array([0.197, 0.803])
    assert phase_fractions == pytest.approx(
        F_expected,
        rel=2e-2
    )

    ch4_gas_fraction = state.speciesAmount('CH4(g)') / phase_fractions[1]
    co2_gas_fraction = state.speciesAmount('CO2(g)') / phase_fractions[1]
    equilibrium_gas_composition = np.array([
        ch4_gas_fraction,
        co2_gas_fraction
    ])
    expected_gas_composition = np.array([0.488, 0.512])
    assert equilibrium_gas_composition == pytest.approx(
        expected_gas_composition,
        rel=1e-2  # 1% rel error
    )

    ch4_liq_fraction = state.speciesAmount('CH4(liq)') / phase_fractions[0]
    co2_liq_fraction = state.speciesAmount('CO2(liq)') / phase_fractions[0]
    equilibrium_liq_composition = np.array([
        ch4_liq_fraction,
        co2_liq_fraction
    ])
    expected_liq_composition = np.array([0.042, 0.958])
    assert equilibrium_liq_composition == pytest.approx(
        expected_liq_composition,
        rel=1e-2  # 1% rel error
    )
