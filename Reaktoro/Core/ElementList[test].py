# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
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
# along with this library. If not, see <http:#www.gnu.org/licenses/>.


from reaktoro import *
import pytest


def testElementList():
    #-------------------------------------------------------------------------
    # TESTING CONSTRUCTOR: ElementList(formulas)
    #-------------------------------------------------------------------------
    elements = ElementList([
        Element("H" , 0.001007940, "Hydrogen" , ["tag1"]),
        Element("He", 0.004002602, "Helium"   , ["tag1"]),
        Element("Li", 0.006941000, "Lithium"  , ["tag1"]),
        Element("Be", 0.009012180, "Beryllium", ["tag2"]),
        Element("B" , 0.010811000, "Boron"    , ["tag2"]),
        Element("C" , 0.012011000, "Carbon"   , ["tag2"]),
        Element("N" , 0.014006740, "Nitrogen" , ["tag3"]),
        Element("O" , 0.015999400, "Oxygen"   , ["tag3"]),
        Element("F" , 0.018998403, "Fluorine" , ["tag1", "tag3"]),
    ])

    assert elements.size() == 9

    assert elements[0].symbol()  == "H"
    assert elements[1].symbol()  == "He"
    assert elements[2].symbol()  == "Li"
    assert elements[3].symbol()  == "Be"
    assert elements[4].symbol()  == "B"
    assert elements[5].symbol()  == "C"
    assert elements[6].symbol()  == "N"
    assert elements[7].symbol()  == "O"
    assert elements[8].symbol()  == "F"

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::find
    #-------------------------------------------------------------------------
    assert elements.find("H")  == 0
    assert elements.find("He") == 1
    assert elements.find("Li") == 2
    assert elements.find("Be") == 3
    assert elements.find("B")  == 4
    assert elements.find("C")  == 5
    assert elements.find("N")  == 6
    assert elements.find("O")  == 7
    assert elements.find("F")  == 8

    assert elements.find("Xy") >= elements.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::findWithSymbol
    #-------------------------------------------------------------------------
    assert elements.findWithSymbol("H")  == 0
    assert elements.findWithSymbol("He") == 1
    assert elements.findWithSymbol("Li") == 2
    assert elements.findWithSymbol("Be") == 3
    assert elements.findWithSymbol("B")  == 4
    assert elements.findWithSymbol("C")  == 5
    assert elements.findWithSymbol("N")  == 6
    assert elements.findWithSymbol("O")  == 7
    assert elements.findWithSymbol("F")  == 8

    assert elements.findWithSymbol("Xy") >= elements.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::findWithName
    #-------------------------------------------------------------------------
    assert elements.findWithName("Hydrogen")  == 0
    assert elements.findWithName("Helium")    == 1
    assert elements.findWithName("Lithium")   == 2
    assert elements.findWithName("Beryllium") == 3
    assert elements.findWithName("Boron")     == 4
    assert elements.findWithName("Carbon")    == 5
    assert elements.findWithName("Nitrogen")  == 6
    assert elements.findWithName("Oxygen")    == 7
    assert elements.findWithName("Fluorine")  == 8

    assert elements.findWithName("Xyrium") >= elements.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::index
    #-------------------------------------------------------------------------
    assert elements.index("H")  == 0
    assert elements.index("He") == 1
    assert elements.index("Li") == 2
    assert elements.index("Be") == 3
    assert elements.index("B")  == 4
    assert elements.index("C")  == 5
    assert elements.index("N")  == 6
    assert elements.index("O")  == 7
    assert elements.index("F")  == 8

    with pytest.raises(Exception):
        elements.index("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::indexWithSymbol
    #-------------------------------------------------------------------------
    assert elements.indexWithSymbol("H")  == 0
    assert elements.indexWithSymbol("He") == 1
    assert elements.indexWithSymbol("Li") == 2
    assert elements.indexWithSymbol("Be") == 3
    assert elements.indexWithSymbol("B")  == 4
    assert elements.indexWithSymbol("C")  == 5
    assert elements.indexWithSymbol("N")  == 6
    assert elements.indexWithSymbol("O")  == 7
    assert elements.indexWithSymbol("F")  == 8

    with pytest.raises(Exception):
        elements.indexWithSymbol("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::indexWithName
    #-------------------------------------------------------------------------
    assert elements.indexWithName("Hydrogen")  == 0
    assert elements.indexWithName("Helium")    == 1
    assert elements.indexWithName("Lithium")   == 2
    assert elements.indexWithName("Beryllium") == 3
    assert elements.indexWithName("Boron")     == 4
    assert elements.indexWithName("Carbon")    == 5
    assert elements.indexWithName("Nitrogen")  == 6
    assert elements.indexWithName("Oxygen")    == 7
    assert elements.indexWithName("Fluorine")  == 8

    with pytest.raises(Exception):
        elements.indexWithName("Xyrium")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::get
    #-------------------------------------------------------------------------
    assert elements.get("H").symbol()  == "H"
    assert elements.get("He").symbol() == "He"
    assert elements.get("Li").symbol() == "Li"
    assert elements.get("Be").symbol() == "Be"
    assert elements.get("B").symbol()  == "B"
    assert elements.get("C").symbol()  == "C"
    assert elements.get("N").symbol()  == "N"
    assert elements.get("O").symbol()  == "O"
    assert elements.get("F").symbol()  == "F"

    with pytest.raises(Exception):
        elements.get("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::getWithSymbol
    #-------------------------------------------------------------------------
    assert elements.getWithSymbol("H") .symbol() == "H"
    assert elements.getWithSymbol("He").symbol() == "He"
    assert elements.getWithSymbol("Li").symbol() == "Li"
    assert elements.getWithSymbol("Be").symbol() == "Be"
    assert elements.getWithSymbol("B") .symbol() == "B"
    assert elements.getWithSymbol("C") .symbol() == "C"
    assert elements.getWithSymbol("N") .symbol() == "N"
    assert elements.getWithSymbol("O") .symbol() == "O"
    assert elements.getWithSymbol("F") .symbol() == "F"

    with pytest.raises(Exception):
        elements.getWithSymbol("Xy")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::getWithName
    #-------------------------------------------------------------------------
    assert elements.getWithName("Hydrogen") .name() == "Hydrogen"
    assert elements.getWithName("Helium")   .name() == "Helium"
    assert elements.getWithName("Lithium")  .name() == "Lithium"
    assert elements.getWithName("Beryllium").name() == "Beryllium"
    assert elements.getWithName("Boron")    .name() == "Boron"
    assert elements.getWithName("Carbon")   .name() == "Carbon"
    assert elements.getWithName("Nitrogen") .name() == "Nitrogen"
    assert elements.getWithName("Oxygen")   .name() == "Oxygen"
    assert elements.getWithName("Fluorine") .name() == "Fluorine"

    with pytest.raises(Exception):
        elements.getWithName("Xyrium")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withSymbols
    #-------------------------------------------------------------------------
    filtered = elements.withSymbols("He O C")

    assert filtered.size() == 3
    assert filtered.indexWithSymbol("He") < filtered.size()
    assert filtered.indexWithSymbol("O" ) < filtered.size()
    assert filtered.indexWithSymbol("C" ) < filtered.size()

    with pytest.raises(Exception):
        elements.withSymbols("He ABC O XYZ C")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withNames
    #-------------------------------------------------------------------------
    filtered = elements.withNames("Helium Oxygen Carbon")

    assert filtered.size() == 3
    assert filtered.indexWithName("Helium") < filtered.size()
    assert filtered.indexWithName("Oxygen") < filtered.size()
    assert filtered.indexWithName("Carbon") < filtered.size()

    with pytest.raises(Exception):
        elements.withSymbols("Helium Xylium")

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withTag
    #-------------------------------------------------------------------------
    filtered = elements.withTag("tag1")

    assert filtered.size() == 4
    assert filtered.indexWithSymbol("H" ) < filtered.size()
    assert filtered.indexWithSymbol("He") < filtered.size()
    assert filtered.indexWithSymbol("Li") < filtered.size()
    assert filtered.indexWithSymbol("F" ) < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withoutTag
    #-------------------------------------------------------------------------
    filtered = elements.withoutTag("tag1")

    assert filtered.size() == 5
    assert filtered.indexWithSymbol("Be") < filtered.size()
    assert filtered.indexWithSymbol("B" ) < filtered.size()
    assert filtered.indexWithSymbol("C" ) < filtered.size()
    assert filtered.indexWithSymbol("N" ) < filtered.size()
    assert filtered.indexWithSymbol("O" ) < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withTags
    #-------------------------------------------------------------------------
    filtered = elements.withTags(["tag1", "tag2"])
    assert filtered.size() == 0

    filtered = elements.withTags(["tag1", "tag3"])
    assert filtered.size() == 1
    assert filtered.index("F") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::withoutTags
    #-------------------------------------------------------------------------
    filtered = elements.withTags(["tag1", "tag2", "tag3"])

    assert filtered.size() == 0

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::append
    #-------------------------------------------------------------------------
    elements.append(Element().withSymbol("Xy").withName("Xyrium"))

    assert elements.indexWithSymbol("Xy")   < elements.size()
    assert elements.indexWithName("Xyrium") < elements.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: ElementList::__iter__
    #-------------------------------------------------------------------------
    for i, element in enumerate(elements):
        assert element.name() == elements[i].name()
