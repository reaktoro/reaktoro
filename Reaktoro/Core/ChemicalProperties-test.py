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
from pytest import approx

from reaktoro import ChemicalEditor, ChemicalSystem


@pytest.fixture
def chemical_system():
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--".split())
    editor.addGaseousPhase("H2O(g) CO2(g)".split())
    editor.addMineralPhase("Graphite")
    return ChemicalSystem(editor)


@pytest.fixture
def chemical_properties(chemical_system):
    # A sensible value for temperature (in K)
    T = 300

    # A sensible value for pressure (in Pa)
    P = 1e5

    # A sensible array of species amounts
    n = np.array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0])

    return chemical_system.properties(T, P, n)


def test_chemical_properties_subvolume(chemical_properties):
    phases = [0]
    volume = chemical_properties.subvolume(phases)
    assert volume.val == approx(0.0010136929)
