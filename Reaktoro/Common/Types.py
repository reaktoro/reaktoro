# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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
import pytest


def testTypes():

    indices = Indices()

    assert indices.size() == 0

    indices.append(1)

    assert indices.size() == 1
    assert len(indices) == 1

    indices.push_back(2)

    assert indices.size() == 2
    assert len(indices) == 2

    assert indices.front() == 1
    assert indices.back() == 2

    indices.append(3)

    assert indices[0] == 1
    assert indices[1] == 2
    assert indices[2] == 3

    for i, val in enumerate(indices):
        assert val == i + 1

