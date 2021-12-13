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


def test_euniquac_params_DTU_values_initialization(ri_values_dtu):
    """
    Test if ri's DTU default parameters are properly initialized and stored. Moreover, the getter
    EUNIQUACParams.ri() is also tested.
    """
    euniquac_params = rkt.EUNIQUACParams()
    euniquac_params.setDTUvalues()
    ri_values = euniquac_params.ri()
    assert type(ri_values) == dict

    ri_expected = ri_values_dtu
    for species in ri_expected:
        assert ri_values[species] == ri_expected[species]


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
