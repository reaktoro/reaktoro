# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pytest

import reaktoro as rkt


@pytest.fixture
def species_names():
    species = [
        "H2O(l)",
        "H+",
        "Na+",
        "K+",
        "NH4+",
        "Cl-",
        "SO4--",
        "HSO4-",
        "NO3-",
        "OH-",
        "CO3--",
        "HCO3-",
        "NN-",
        "NN+",
        "Ca++",
        "Sr++",
        "Ba++",
    ]
    return species


@pytest.fixture
def ri_values_dtu(species_names):
    ri = [
        0.92,
        0.13779,
        1.4034,
        2.2304,
        4.8154,
        10.386,
        12.794,
        19.588,
        5.4041,
        9.3973,
        11.09,
        3.755,
        0.13779,
        0.13779,
        3.87,
        4.1852,
        1.525,
    ]
    ri_values = {}
    for i, species in enumerate(species_names):
        ri_values[species] = ri[i]

    return ri_values


@pytest.fixture
def qi_values_dtu(species_names):
    qi = [
        1.4,
        1.00e-15,
        1.199,
        2.4306,
        4.6028,
        10.197,
        12.444,
        22.495,
        6.2074,
        8.8171,
        11.32,
        4.712,
        1.00e-15,
        1.00e-15,
        1.48,
        0.2818,
        0.2278,
    ]
    qi_values = {}
    for i, species in enumerate(species_names):
        qi_values[species] = qi[i]

    return qi_values


def test_species_id_map_setup():
    """
    Test if the species id map can be properly set and retrieved.
    """
    bips_species_id_map = {
        "Na+": 0,
        "Cl-": 1,
    }

    uniquac_params = rkt.EUNIQUACParams()
    uniquac_params.bipsSpeciesIds(bips_species_id_map)
    assert uniquac_params.bipsSpeciesIds() == bips_species_id_map


def test_euniquac_params_DTU_values_initialization(ri_values_dtu, qi_values_dtu):
    """
    Test if ri and qi DTU default parameters are properly initialized and stored.
    Moreover, the getters EUNIQUACParams.ri(), EUNIQUACParams.qi(), EUNIQUACParams.uij_0(),
    and EUNIQUACParams.uij_T() are also tested.
    """
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.setDTUvalues()
    ri_values = euniquac_params.ri()
    qi_values = euniquac_params.qi()
    assert type(ri_values) == dict
    assert type(qi_values) == dict

    ri_expected = ri_values_dtu
    qi_expected = qi_values_dtu
    for species in ri_expected:
        assert ri_values[species] == ri_expected[species]
        assert qi_values[species] == qi_expected[species]

    bips_species_id_map_dtu = {
        "H2O(l)": 0,
        "H+": 1,
        "Na+": 2,
        "K+": 3,
        "NH4+": 4,
        "Cl-": 5,
        "SO4--": 6,
        "HSO4-": 7,
        "NO3-": 8,
        "OH-": 9,
        "CO3--": 10,
        "HCO3-": 11,
    }
    species_id_stored = euniquac_params.bipsSpeciesIds()
    assert species_id_stored == bips_species_id_map_dtu

    uij_0_dtu = np.array([
        [0,      1e5,     733.286, 535.023, 54.0297, 1523.39, 752.879, 602.252, 998.92, 600.495, 328.141,  118.702],
        [1e5,     0,       1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10, 1.0E+10],
        [733.286, 1.0E+10, 0,       -46.194, 375.977, 1443.23, 845.135, 469.488, 797.474, 1398.14, 476.956, 980.982],
        [535.023, 1.0E+10, -46.194, 0,       184.288, 1465.18, 913.824, 445.673, 818.568, 1805.75, 1370.57, 1123.44],
        [54.0297, 1.0E+10, 375.977, 184.288, 0,       1385.31, 677.178, 418.886, 807.246, 2500,    0,       0],
        [1523.39, 1.0E+10, 1443.23, 1465.18, 1385.31, 2214.81, 2036.06, 0,       2175.02, 1895.52, 2372.94, 2014.18],
        [752.879, 1.0E+10, 845.135, 913.824, 677.178, 2036.06, 1265.83, 1004.16, 1757.79, 1225.67, 1158.91, 956.609],
        [602.252, 1.0E+10, 469.488, 445.673, 418.886, 0,       1004.16, 822.989, 0,       2500,    2500,     2500],
        [998.92,  1.0E+10, 797.474, 818.568, 807.246, 2175.02, 1757.79, 0,       2753.71, 1379.95, 0,       0,],
        [600.495, 1.0E+10, 1398.14, 1805.75, 2500,    1895.52, 1225.67, 2500,    1379.95, 1562.88, 1339.04, 2500],
        [328.141, 1.0E+10, 476.956, 1370.57, 0,       2372.94, 1158.91, 2500,    0,       1339.04, 1065.97, 565.786],
        [118.702, 1.0E+10, 980.982, 1123.44, 0,       2014.18, 956.609, 2500,    0,       2500,    565.786, 253.461],
    ])
    uij_0_stored = euniquac_params.uij_0()
    assert np.allclose(uij_0_stored, uij_0_dtu)

    uij_T_dtu = np.array([
        [0, 0, 0.4872, 0.9936, 0.5855, 14.631, 9.4905, 5.5499, 9.3251, 8.5455, -0.5059, 0.96],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0.4872, 0, 0, 0.119, -0.292, 15.635, 11.681, 6.5329, 9.3047, 20.278, 2.8191, 13.296],
        [0.9936, 0, 0.119, 0, 1.0985, 15.329, 12.278, 5.2231, 10.302, 27.283, 5.0517, 2.7217],
        [0.5855, 0, -0.292, 1.0985, 0, 14.848, 10.356, 5.1053, 9.9705, 0, 0, 0],
        [14.631, 0, 15.635, 15.329, 14.848, 14.436, 12.407, 0, 13.449, 13.628, 8.9402, 11.636],
        [9.4905, 0, 11.681, 12.278, 10.356, 12.407, 8.3194, 8.7877, 6.9013, 8.5902, 6.7255, 8.6656],
        [5.5499, 0, 6.5329, 5.2231, 5.1053, 0, 8.7877, 5.7939, 0, 0, 0, 0],
        [9.3251, 0, 9.3047, 10.302, 9.9705, 13.449, 6.9013, 0, 2.2866, 6.6369, 0, 0],
        [8.5455, 0, 20.278, 27.283, 0, 13.628, 8.5902, 0, 6.6369, 5.6169, 3.8129, 0],
        [-0.5059, 0, 2.8191, 5.0517, 0, 8.9402, 6.7255, 0, 0, 3.8129, -4.4653, -1.2033],
        [0.96, 0, 13.296, 2.7217, 0, 11.636, 8.6656, 0, 0, 0, -1.2033, 1.2824]
    ])
    uij_T_stored = euniquac_params.uij_T()
    assert np.allclose(uij_T_stored, uij_T_dtu)


def test_euniquac_ri_setup():
    """
    Test if the E-UNIQUAC ri's setters and getters are working properly.
    """
    ri_ion_1 = 300.0
    new_ion_1_dict = {"ion1": ri_ion_1}
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.ri(new_ion_1_dict)
    assert euniquac_params.ri("ion1") == new_ion_1_dict["ion1"]

    ri_ion_2_name = "ion2"
    ri_ion_2_value = 400.0
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.ri(ri_ion_2_name, ri_ion_2_value)
    assert euniquac_params.ri(ri_ion_2_name) == ri_ion_2_value


def test_euniquac_qi_setup():
    """
    Test if the E-UNIQUAC qi's setters and getters are working properly.
    """
    qi_ion_1 = 300.0
    new_ion_1_dict = {"ion1": qi_ion_1}
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.qi(new_ion_1_dict)
    assert euniquac_params.qi("ion1") == new_ion_1_dict["ion1"]

    qi_ion_2_name = "ion2"
    qi_ion_2_value = 400.0
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.ri(qi_ion_2_name, qi_ion_2_value)
    assert euniquac_params.ri(qi_ion_2_name) == qi_ion_2_value


def test_uij_bips_setup():
    """
    Test if energetic BIPs matrices are properly working for getters and setters.
    """
    # Please note that the BIPs matrices should be symmetric
    uij_0_custom = np.array([
        [1.0, 2.0],
        [2.0, 3.0],
    ])
    uij_T_custom = np.array([
        [4.0, 5.0],
        [5.0, 6.0],
    ])
    bips_species_id_map = {
        "Na+": 0,
        "Cl-": 1,
    }

    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.set_uij_bips(uij_0_custom, uij_T_custom, bips_species_id_map)

    uij_0_stored = euniquac_params.uij_0()
    uij_T_stored = euniquac_params.uij_T()
    assert euniquac_params.bipsSpeciesIds() == bips_species_id_map
    assert np.allclose(uij_0_stored, uij_0_custom)
    assert np.allclose(uij_T_stored, uij_T_custom)
    assert euniquac_params.uij_0("Na+", "Na+") == 1.0
    assert euniquac_params.uij_0("Na+", "Cl-") == 2.0
    assert euniquac_params.uij_0("Cl-", "Cl-") == 3.0
    assert euniquac_params.uij_T("Na+", "Na+") == 4.0
    assert euniquac_params.uij_T("Na+", "Cl-") == 5.0
    assert euniquac_params.uij_T("Cl-", "Cl-") == 6.0


def test_bips_values_update():
    """
    Test if energy BIPs can be modified after a first setup.
    """
    uij_0_custom = np.array([
        [1.0, 2.0],
        [2.0, 3.0],
    ])
    uij_T_custom = np.array([
        [4.0, 5.0],
        [5.0, 6.0],
    ])
    bips_species_id_map = {
        "Na+": 0,
        "Cl-": 1,
    }

    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.set_uij_bips(uij_0_custom, uij_T_custom, bips_species_id_map)

    euniquac_params.uij_0("Na+", "Cl-", 3.5)
    assert euniquac_params.uij_0("Na+", "Cl-") == 3.5
    assert euniquac_params.uij_0("Cl-", "Na+") == 3.5

    euniquac_params.uij_T("Na+", "Cl-", 1.5)
    assert euniquac_params.uij_T("Na+", "Cl-") == 1.5
    assert euniquac_params.uij_T("Cl-", "Na+") == 1.5

    euniquac_params.uij_0("Na+", "Na+", 2.5)
    euniquac_params.uij_T("Na+", "Na+", 2.5)
    assert euniquac_params.uij_0("Na+", "Na+") == 2.5
    assert euniquac_params.uij_T("Na+", "Na+") == 2.5


def test_add_new_species_params():
    uij_0_custom = np.array([
        [1.0, 2.0],
        [2.0, 3.0],
    ])
    uij_T_custom = np.array([
        [4.0, 5.0],
        [5.0, 6.0],
    ])
    bips_species_id_map = {
        "Na+": 0,
        "Cl-": 1,
    }

    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.set_uij_bips(uij_0_custom, uij_T_custom, bips_species_id_map)

    # New "dummy" addition. Be careful to add a species that can be handled internally by Reaktoro.
    new_species_name = "K+"
    q_value = 2.0
    r_value = 1.0
    new_uij_0_values = {"Na+": 3.5}
    new_uij_T_values = {"Cl-": -2.5}
    euniquac_params.addNewSpeciesParameters(
        new_species_name,
        q_value,
        r_value,
        new_uij_0_values,
        new_uij_T_values
    )

    # Check if all previous values are untouched
    assert euniquac_params.uij_0("Na+", "Na+") == 1.0
    assert euniquac_params.uij_0("Na+", "Cl-") == 2.0
    assert euniquac_params.uij_0("Cl-", "Cl-") == 3.0
    assert euniquac_params.uij_T("Na+", "Na+") == 4.0
    assert euniquac_params.uij_T("Na+", "Cl-") == 5.0
    assert euniquac_params.uij_T("Cl-", "Cl-") == 6.0

    # Check the new values between the new species and the ones already defined (for uij_0)
    default_u0_value = 2500.0
    assert euniquac_params.uij_0(new_species_name, new_species_name) == 0.0
    assert euniquac_params.uij_0("Na+", new_species_name) == 3.5
    assert euniquac_params.uij_0(new_species_name, "Na+") == 3.5
    assert euniquac_params.uij_0("Cl-", new_species_name) == default_u0_value
    assert euniquac_params.uij_0(new_species_name, "Cl-") == default_u0_value

    # Check the new values between the new species and the ones already defined (for uij_T)
    default_uT_value = 0.0
    assert euniquac_params.uij_T("Cl-", new_species_name) == -2.5
    assert euniquac_params.uij_T(new_species_name, "Cl-") == -2.5
    assert euniquac_params.uij_T("Na+", new_species_name) == default_uT_value
    assert euniquac_params.uij_T(new_species_name, "Na+") == default_uT_value
    assert euniquac_params.uij_T(new_species_name, new_species_name) == default_uT_value


def test_euniquac_fallback_lr_only_flag():
    """
    Test if the E-UNIQUAC fallback to long range only is working properly
    """

    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.setLongRangeOnlyForSpeciesMissingParameters()
    assert True


def test_euniquac_fallback_get_species_missing_euniquac_params():
    """
    Test if the E-UNIQUAC fallback to long range only is working properly
    """

    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.setLongRangeOnlyForSpeciesMissingParameters()
    assert True

def test_euniquac_fallback_all_species_are_known_match_previous_results():
    """
    Test if the E-UNIQUAC with fallback modifications but with all species
    with known q, r and bips match the previous results
    """

    list_aqueous_species = [
        "H2O(l)",
        "H+",
        "OH-",
        "Na+",
        'Cl-',
    ]
    gas_species = []

    mineral_name = 'Halite'

    # Create Euniquac
    editor = rkt.ChemicalEditor()
    editor.addMineralPhase(mineral_name)
    aqueous_phase = editor.addAqueousPhase(list_aqueous_species)

    if len(gas_species):
        gas_phase = editor.addGaseousPhase(gas_species)
        gas_phase.setChemicalModelSoaveRedlichKwong()

    euniquac_params = rkt.EUNIQUACParams()
    editor.aqueousPhase().setChemicalModelEUNIQUAC(euniquac_params)
    system_euniquac = rkt.ChemicalSystem(editor)


    problem = rkt.EquilibriumProblem(system_euniquac)
    problem.setTemperature(25.0, "celsius")
    problem.setPressure(1.0, "atm")
    
    problem.add("H2O", 1, "kg")
    problem.add(mineral_name, 100, "mol") # excess quantity

    state = rkt.equilibrate(problem)

    solubility = state.elementAmountInPhase('Na', 'Aqueous')*1e3
    pH = rkt.ChemicalProperty.pH(system_euniquac)(state.properties()).val
    mol_species = state.speciesAmounts

    assert np.isclose(solubility, 7305.134115874772)
    assert np.isclose(pH, 8.001947888594076)
    expected_mols = [6.55084350e+01, 6.14295491e-08, 6.14295491e-08, 7.30513412e+00, 7.30513412e+00, 9.26948659e+01]
    assert np.testing.assert_array_almost_equal(mol_species, expected_mols)