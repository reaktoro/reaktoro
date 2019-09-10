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

import numpy as np
from reaktoro import ChemicalProperties


def test_chemical_properties_subvolume(chemical_properties):
    phases = [0]
    volume = chemical_properties.subvolume(phases)
    assert volume.val == pytest.approx(0.0010136929)


def test_chemical_properties_not_updated(chemical_system):
    chemical_properties = ChemicalProperties(chemical_system)
    assert np.isnan(chemical_properties.temperature().val)
    assert np.isnan(chemical_properties.pressure().val)

    with pytest.raises(RuntimeError, match=r"Cannot proceed with method ChemicalProperties::update."):
        chemical_properties.update(np.array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0]))

