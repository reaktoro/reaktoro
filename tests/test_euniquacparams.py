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

import reaktoro as rkt


def test_euniquac_ri_setup():
    """
    Test if the E-UNIQUAC ri's setters and getters are working properly.
    """
    editor = rkt.ChemicalEditor()
    editor.addGaseousPhase(["H2O(g)"])
    editor.addAqueousPhase(["H2O(l)", "H+", "OH-", "Na+", "Cl-"])

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
