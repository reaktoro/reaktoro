import numpy as np

from reaktoro import (
    ChemicalEditor,
    ChemicalState,
    ChemicalSystem,
    Database,
    EquilibriumSolver,
    BinaryInteractionParams,
    CubicEOSParams
)


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
    assert np.allclose(bips.k, k)
    assert np.allclose(bips.kT, kT)
    assert np.allclose(bips.kTT, kTT)


def test_bips_calculation_function():
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
    
    assert np.allclose(bips_calculated.k, bips_expected.k)
    assert np.allclose(bips_calculated.kT, bips_expected.kT)
    assert np.allclose(bips_calculated.kTT, bips_expected.kTT)


def test_bips_setup_without_derivatives():
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
    
    assert np.allclose(bips_calculated.k, bips_expected.k)
