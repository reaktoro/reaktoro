# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
#
# This library is free software you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *
import pytest


def testCriticalProps():

    assert CriticalProps.size() == len(CriticalProps.data())

    assert CriticalProps.find("CO2") < CriticalProps.size()
    assert CriticalProps.find("H2O") < CriticalProps.size()

    assert CriticalProps.find("XYZ") == CriticalProps.size()

    crprops = SubstanceCriticalProps("XYZ")
    crprops.setTemperature(123.4, "K")
    crprops.setPressure(123.0, "bar")
    crprops.setAcentricFactor(0.1234)

    CriticalProps.append(crprops)

    assert CriticalProps.find("XYZ") < CriticalProps.size()

    assert CriticalProps.get("CO2").temperature()    == 304.20
    assert CriticalProps.get("CO2").pressure()       == 73.83e+5
    assert CriticalProps.get("CO2").acentricFactor() == 0.2240

    CriticalProps.setMissingAs("He")

    assert CriticalProps.get("ABCDEF").temperature()    == CriticalProps.get("He").temperature()
    assert CriticalProps.get("ABCDEF").pressure()       == CriticalProps.get("He").pressure()
    assert CriticalProps.get("ABCDEF").acentricFactor() == CriticalProps.get("He").acentricFactor()

    CriticalProps.setMissingAs("H2")

    assert CriticalProps.get("ABCDEF").temperature()    == CriticalProps.get("H2").temperature()
    assert CriticalProps.get("ABCDEF").pressure()       == CriticalProps.get("H2").pressure()
    assert CriticalProps.get("ABCDEF").acentricFactor() == CriticalProps.get("H2").acentricFactor()

    crprops = SubstanceCriticalProps("CO2")
    crprops.setTemperature(123.4, "K")
    crprops.setPressure(123.0, "bar")
    crprops.setAcentricFactor(0.1234)

    CriticalProps.overwrite(crprops)

    assert CriticalProps.get("CO2").temperature()    == 123.4
    assert CriticalProps.get("CO2").pressure()       == 123.0e5
    assert CriticalProps.get("CO2").acentricFactor() == 0.1234
