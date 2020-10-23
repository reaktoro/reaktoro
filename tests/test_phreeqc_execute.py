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
from pytest import approx

def test_phreeqc_execute(shared_datadir):
    """Test function the execute() method from class Phreeqc."""

    # all '\' (from Windows path format) needs to be change to '/'.
    database_path = '/'.join([i.replace('\\', '') for i in (shared_datadir / 'phreeqc.dat').parts])
    input_path = '/'.join([i.replace('\\', '') for i in (shared_datadir / 'IW.pqi').parts])

    p = Phreeqc(database_path)
    p.execute(input_path)

    assert approx(p.temperature()) == 298.15
    assert approx(p.pressure()) == 101324.99997084374
    assert approx(p.numElements()) == 11
    assert approx(p.numSpecies()) == 44
    assert approx(p.numPhases()) == 1

