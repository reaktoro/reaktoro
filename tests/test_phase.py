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

from reaktoro import *


def test_phase():
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- CO2(aq) CO3-- HCO3-")

    system = ChemicalSystem(editor)

    # Test Phase::elements()
    phase = system.phase(0)

    assert len(phase.elements()) == 4

    symbols = [element.name() for element in phase.elements()]

    assert "H" in symbols
    assert "O" in symbols
    assert "C" in symbols
    assert "Z" in symbols
