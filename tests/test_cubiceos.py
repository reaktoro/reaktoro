from numpy.testing import assert_array_equal

from reaktoro import (
    ChemicalEditor,
    ChemicalState,
    ChemicalSystem,
    Database,
    EquilibriumSolver,
    BinaryInteractionParams,
    CubicEOSParams,
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
    k_00 = k_11 = k_22 = 0.0
    k_01 = k_10 = 0.01
    k_02 = k_20 = 0.50
    k_12 = k_21 = 0.40
    k = [
        [k_00, k_01, k_02],
        [k_10, k_11, k_12],
        [k_20, k_21, k_22]
    ]

    kT_00 = kT_11 = kT_22 = 0.0
    kT_01 = kT_10 = 0.0
    kT_02 = kT_20 = 0.0
    kT_12 = kT_21 = 0.0
    kT = [
        [kT_00, kT_01, kT_02],
        [kT_10, kT_11, kT_12],
        [kT_20, kT_21, kT_22]
    ]

    kTT_00 = kTT_11 = kTT_22 = 0.0
    kTT_01 = kTT_10 = 0.0
    kTT_02 = kTT_20 = 0.0
    kTT_12 = kTT_21 = 0.0
    kTT = [
        [kTT_00, kTT_01, kTT_02],
        [kTT_10, kTT_11, kTT_12],
        [kTT_20, kTT_21, kTT_22]
    ]

    bips = BinaryInteractionParams(k, kT, kTT)

    assert_array_equal(bips.k, k)
    assert_array_equal(bips.kT, kT)
    assert_array_equal(bips.kTT, kTT)


def test_bips_calculation_function():
    """
    Test the wrapper for the function that calculates the BIPs for mixtures
    modeled by a Cubic EOS. In this particular case, the BIPs depend on the
    temperature.
    """
    def bips_function(T):
        k_00 = k_11 = k_22 = 0.0 * T
        k_01 = k_10 = 0.01 * T
        k_02 = k_20 = 0.50 * T
        k_12 = k_21 = 0.40 * T
        k = [
            [k_00, k_01, k_02],
            [k_10, k_11, k_12],
            [k_20, k_21, k_22]
        ]

        kT_00 = kT_11 = kT_22 = 0.0
        kT_01 = kT_10 = 0.01
        kT_02 = kT_20 = 0.50
        kT_12 = kT_21 = 0.40
        kT = [
            [kT_00, kT_01, kT_02],
            [kT_10, kT_11, kT_12],
            [kT_20, kT_21, kT_22]
        ]

        kTT_00 = kTT_11 = kTT_22 = 0.0
        kTT_01 = kTT_10 = 0.0
        kTT_02 = kTT_20 = 0.0
        kTT_12 = kTT_21 = 0.0
        kTT = [
            [kTT_00, kTT_01, kTT_02],
            [kTT_10, kTT_11, kTT_12],
            [kTT_20, kTT_21, kTT_22]
        ]
        bips_values = BinaryInteractionParams(k, kT, kTT)
        return bips_values

    T_dummy = 2.0
    bips_expected = bips_function(T_dummy)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T_dummy)

    assert_array_equal(bips_calculated.k, bips_expected.k)
    assert_array_equal(bips_calculated.kT, bips_expected.kT)
    assert_array_equal(bips_calculated.kTT, bips_expected.kTT)


def test_bips_setup_without_derivatives():
    """
    This test takes into account the case where the BIPs are constants. In such a case, define
    matrices for kT and kTT are not necessary since all entries would be zero. Thus, only k
    should be defined.
    """
    def bips_function(T):
        k_00 = k_11 = k_22 = 0.0
        k_01 = k_10 = 0.01
        k_02 = k_20 = 0.50
        k_12 = k_21 = 0.40
        k = [
            [k_00, k_01, k_02],
            [k_10, k_11, k_12],
            [k_20, k_21, k_22]
        ]
        bips_values = BinaryInteractionParams(k)
        return bips_values

    T_dummy = 2.0
    bips_expected = bips_function(T_dummy)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T_dummy)

    assert_array_equal(bips_calculated.k, bips_expected.k)    
    assert len(bips_calculated.kT) == 0
    assert len(bips_calculated.kTT) == 0
