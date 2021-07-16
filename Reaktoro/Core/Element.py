# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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


def testElement():
    # When using constructor Element()
    element = Element()

    element = element.withSymbol("Na")
    assert element.symbol() == "Na"

    element = element.withName("Sodium")
    assert element.name() == "Sodium"

    element = element.withMolarMass(0.022989768)
    assert element.molarMass() == 0.022989768

    element = element.withTags(["tag1", "tag2"])
    assert element.tags() == ["tag1", "tag2"]

    # When using constructor Element(Element::Attribs)
    element = Element(symbol="H", molar_mass=0.001007940)

    assert element.symbol() == "H"
    assert element.molarMass() == 0.001007940
    assert element.name() == "H"
    assert element.tags() == []

    element = Element(symbol="H", molar_mass=0.001007940, name="Hydrogen", tags=["tag1", "tag2", "tag3"])
    assert element.name() == "Hydrogen"
    assert element.tags() == ["tag1", "tag2", "tag3"]

    with pytest.raises(Exception):
        element.withMolarMass(-1.0)

