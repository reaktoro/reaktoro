import numpy as np
from numpy.testing import assert_array_equal

from reaktoro import (
    ChemicalEditor,
    ChemicalState,
    ChemicalSystem,
    Database,
    EquilibriumSolver,
    BinaryInteractionParams,
    CubicEOSParams,
    ThermoScalar
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
    k_00 = k_11 = k_22 = ThermoScalar(0.0)
    k_01 = k_10 = ThermoScalar(0.01)
    k_02 = k_20 = ThermoScalar(0.50)
    k_12 = k_21 = ThermoScalar(0.40)
    k = [
        [k_00, k_01, k_02],
        [k_10, k_11, k_12],
        [k_20, k_21, k_22]
    ]

    kT_00 = kT_11 = kT_22 = ThermoScalar(0.0)
    kT_01 = kT_10 = ThermoScalar(0.0)
    kT_02 = kT_20 = ThermoScalar(0.0)
    kT_12 = kT_21 = ThermoScalar(0.0)
    kT = [
        [kT_00, kT_01, kT_02],
        [kT_10, kT_11, kT_12],
        [kT_20, kT_21, kT_22]
    ]

    kTT_00 = kTT_11 = kTT_22 = ThermoScalar(0.0)
    kTT_01 = kTT_10 = ThermoScalar(0.0)
    kTT_02 = kTT_20 = ThermoScalar(0.0)
    kTT_12 = kTT_21 = ThermoScalar(0.0)
    kTT = [
        [kTT_00, kTT_01, kTT_02],
        [kTT_10, kTT_11, kTT_12],
        [kTT_20, kTT_21, kTT_22]
    ]

    bips = BinaryInteractionParams(k, kT, kTT)

    k_expected = np.array([
        [k_00.val, k_01.val, k_02.val],
        [k_10.val, k_11.val, k_12.val],
        [k_20.val, k_21.val, k_22.val]
    ])

    kT_expected = np.array([
        [kT_00.val, kT_01.val, kT_02.val],
        [kT_10.val, kT_11.val, kT_12.val],
        [kT_20.val, kT_21.val, kT_22.val]
    ])

    kTT_expected = np.array([
        [kTT_00.val, kTT_01.val, kTT_02.val],
        [kTT_10.val, kTT_11.val, kTT_12.val],
        [kTT_20.val, kTT_21.val, kTT_22.val]
    ])

    bips_k = np.array(bips.k)
    bips_kT = np.array(bips.kT)
    bips_kTT = np.array(bips.kTT)

    nspecies = bips_k.shape[0]
    for i in range(nspecies):
        for j in range(nspecies):
            assert bips_k[i, j].val == k_expected[i, j]
            assert bips_kT[i, j].val == kT_expected[i, j]
            assert bips_kTT[i, j].val == kTT_expected[i, j]


def test_bips_calculation_function():
    """
    Test the wrapper for the function that calculates the BIPs for mixtures
    modeled by a Cubic EOS. In this particular case, the BIPs depend on the
    temperature.
    """
    def bips_function(T):
        k_00 = k_11 = k_22 = ThermoScalar(0.0 * T.val)
        k_01 = k_10 = ThermoScalar(0.01 * T.val)
        k_02 = k_20 = ThermoScalar(0.50 * T.val)
        k_12 = k_21 = ThermoScalar(0.40 * T.val)
        k = [
            [k_00, k_01, k_02],
            [k_10, k_11, k_12],
            [k_20, k_21, k_22]
        ]

        kT_00 = kT_11 = kT_22 = ThermoScalar(0.0)
        kT_01 = kT_10 = ThermoScalar(0.01)
        kT_02 = kT_20 = ThermoScalar(0.50)
        kT_12 = kT_21 = ThermoScalar(0.40)
        kT = [
            [kT_00, kT_01, kT_02],
            [kT_10, kT_11, kT_12],
            [kT_20, kT_21, kT_22]
        ]

        kTT_00 = kTT_11 = kTT_22 = ThermoScalar(0.0)
        kTT_01 = kTT_10 = ThermoScalar(0.0)
        kTT_02 = kTT_20 = ThermoScalar(0.0)
        kTT_12 = kTT_21 = ThermoScalar(0.0)
        kTT = [
            [kTT_00, kTT_01, kTT_02],
            [kTT_10, kTT_11, kTT_12],
            [kTT_20, kTT_21, kTT_22]
        ]
        bips_values = BinaryInteractionParams(k, kT, kTT)
        return bips_values

    T_dummy = 2.0
    T = ThermoScalar(T_dummy)
    bips_expected = bips_function(T)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T)

    bips_k_calculated = np.array(bips_calculated.k)
    bips_kT_calculated = np.array(bips_calculated.kT)
    bips_kTT_calculated = np.array(bips_calculated.kTT)

    bips_k_expected = np.array(bips_expected.k)
    bips_kT_expected = np.array(bips_expected.kT)
    bips_kTT_expected = np.array(bips_expected.kTT)

    nspecies = bips_k_calculated.shape[0]
    for i in range(nspecies):
        for j in range(nspecies):
            assert bips_k_calculated[i, j].val == bips_k_expected[i, j].val
            assert bips_kT_calculated[i, j].val == bips_kT_expected[i, j].val
            assert bips_kTT_calculated[i, j].val == bips_kTT_expected[i, j].val


def test_bips_setup_without_derivatives():
    """
    This test takes into account the case where the BIPs are constants. In such a case, define
    matrices for kT and kTT are not necessary since all entries would be zero. Thus, only k
    should be defined.
    """
    def bips_function(T):
        k_00 = k_11 = k_22 = ThermoScalar(0.0)
        k_01 = k_10 = ThermoScalar(0.01)
        k_02 = k_20 = ThermoScalar(0.50)
        k_12 = k_21 = ThermoScalar(0.40)
        k = [
            [k_00, k_01, k_02],
            [k_10, k_11, k_12],
            [k_20, k_21, k_22]
        ]
        bips_values = BinaryInteractionParams(k)
        return bips_values

    T_dummy = 2.0
    T = ThermoScalar(T_dummy)
    bips_expected = bips_function(T)

    cubic_eos_params = CubicEOSParams(binary_interaction_values=bips_function)
    bips_calculated = cubic_eos_params.binary_interaction_values(T)

    bips_k_calculated = np.array(bips_calculated.k)

    bips_k_expected = np.array(bips_expected.k)

    nspecies = bips_k_calculated.shape[0]
    for i in range(nspecies):
        for j in range(nspecies):
            assert bips_k_calculated[i, j].val == bips_k_expected[i, j].val
    
    assert len(bips_calculated.kT) == 0
    assert len(bips_calculated.kTT) == 0
