import numpy as np
import pytest

import reaktoro
from reaktoro import (
    ChemicalEditor,
    ChemicalState,
    ChemicalSystem,
    Database,
    EquilibriumSolver,
    BinaryInteractionParams,
    CubicEOSParams
)


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[0] / "data"



def test_CubicEOS_multiple_roots():
    """
    This problem leads to the following CubicEOS roots
    PR - Z1 = 1.00027728
         Z2 = 0.0001655
         Z3 = -0.0011024
    since bmix = 1.635e-05 -> Z3 is an invalid root 
    and since Z3 < Z2 < Z1 -> Z2 is an invalid root.
    Reaktoro should remove Z3, Z2 and proceed instead of removing only Z3 and
    raising the exception "Logic error: it was expected Z roots of size 3, but
    got: 2".
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElementsOf("H2O Fe(OH)2 Fe(OH)3 NH3")
    editor.addGaseousPhaseWithElementsOf("NH3")
    editor.addMineralPhase("Magnetite")

    system = ChemicalSystem(editor)

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    temperature = 298.15
    pressure = 1e5
    b = [
        3.0,
        122.01687012,
        1.0,
        63.50843506,
        0.0,
    ]

    result = solver.approximate(state, temperature, pressure, b)
    assert result.optimum.succeeded

    # check that it doesn't raise an exception
    state.properties()


def test_bips_setup():
    """
    Test the BIPs storage in a BinaryInteractionParams object.
    """
    k = [
        [0.00, 0.01, 0.50],
        [0.01, 0.00, 0.40],
        [0.50, 0.40, 0.00]
    ]
    kT = [
        [0.00, 0.00, 0.00],
        [0.00, 0.00, 0.00],
        [0.00, 0.00, 0.00]
    ]
    kTT = [
        [0.00, 0.00, 0.00],
        [0.00, 0.00, 0.00],
        [0.00, 0.00, 0.00]
    ]

    bips = BinaryInteractionParams(k, kT, kTT)
    assert bips.k.all() == np.array(k).all()
    assert bips.kT.all() == np.array(kT).all()
    assert bips.kTT.all() == np.array(kTT).all()


def test_bips_calculation_function():
    """
    Test the wrapper for the function that calculates the BIPs for mixtures
    modeled by a Cubic EOS. In this particular case, the BIPs depend on the
    temperature.
    """
    def bips_function(T):
        k = [
            [0.00 * T, 0.01 * T, 0.50 * T],
            [0.01 * T, 0.00 * T, 0.40 * T],
            [0.50 * T, 0.40 * T, 0.00 * T]
        ]
        kT = [
            [0.00, 0.01, 0.50],
            [0.01, 0.00, 0.40],
            [0.50, 0.40, 0.00]
        ]
        kTT = [
            [0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00],
            [0.00, 0.00, 0.00]
        ]
        bips_values = BinaryInteractionParams(k, kT, kTT)
        return bips_values

    T_dummy = 2.0
    bips_expected = bips_function(T_dummy)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T_dummy)
    
    assert bips_calculated.k.all() == np.array(bips_expected.k).all()
    assert bips_calculated.kT.all() == np.array(bips_expected.kT).all()
    assert bips_calculated.kTT.all() == np.array(bips_expected.kTT).all()


def test_bips_setup_without_derivatives():
    """
    This test takes into account the case where the BIPs are constants. In such a case, define
    matrices for kT and kTT are not necessary since all entries would be zero. Thus, only k
    should be defined.
    """
    def bips_function(T):
        k = [
            [0.00, 0.01, 0.50],
            [0.01, 0.00, 0.40],
            [0.50, 0.40, 0.00]
        ]
        bips_values = BinaryInteractionParams(k)
        return bips_values

    T_dummy = 2.0
    bips_expected = bips_function(T_dummy)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T_dummy)
    
    assert bips_calculated.k.all() == np.array(bips_expected.k).all()
    assert len(bips_calculated.kT) == 0
    assert len(bips_calculated.kTT) == 0


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
