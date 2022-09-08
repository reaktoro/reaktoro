# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
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


# TODO Implement tests for the python bindings of component Data in Data[test].py
def testData():
    foo = Data()
    foo.add("A", 1.0)
    foo.add("B", 2)

    bar = Data()
    bar.add("C", True)
    bar.add("D", "X")
    bar.add("E", Param(7.0))

    doo = Data()
    doo.add(3.0)
    doo.add(6.0)
    doo.add(9.0)

    params = Data()
    params.add("Foo", foo)
    params.add("Bar", bar)

    assert foo["A"].number() == 1.0
    assert foo["B"].integer() == 2
    assert foo.exists("A") == True
    assert foo.exists("B") == True
    assert foo.exists("C") == False

    assert bar["C"].boolean() == True
    assert bar["D"].string() == "X"
    assert bar["E"].param() == 7.0
    assert bar.exists("C") == True
    assert bar.exists("D") == True
    assert bar.exists("E") == True
    assert bar.exists("F") == False

    assert doo[0].number() == 3.0
    assert doo[1].number() == 6.0
    assert doo[2].number() == 9.0

    assert params["Foo"]["A"].number() == 1.0
    assert params["Foo"]["B"].integer() == 2
    assert params["Bar"]["C"].boolean() == True
    assert params["Bar"]["D"].string() == "X"
    assert params["Bar"]["E"].param() == 7.0
    assert params.exists("Foo") == True
    assert params.exists("Bar") == True
    assert params.exists("Joe") == False
    assert params["Foo"].exists("A") == True
    assert params["Foo"].exists("Z") == False
    assert params["Bar"].exists("C") == True
    assert params["Bar"].exists("Z") == False

