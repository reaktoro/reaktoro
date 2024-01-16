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


def testElementalComposition():

    #-----------------------------------------------------------------
    # Testing constructor ElementalComposition(list[tuple[Element, float]])
    #-----------------------------------------------------------------
    elements = ElementalComposition([
        ( Element("H"), 1.0 ),
        ( Element("C"), 2.0 ),
        ( Element("O"), 3.0 ),
    ])

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.symbols()
    #-----------------------------------------------------------------
    assert elements.symbols() == ["H", "C", "O"]

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.coefficients()
    #-----------------------------------------------------------------
    assert elements.coefficients() == [1.0, 2.0, 3.0]

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.coefficient()
    #-----------------------------------------------------------------
    assert elements.coefficient("H") == 1.0
    assert elements.coefficient("C") == 2.0
    assert elements.coefficient("O") == 3.0
    assert elements.coefficient("X") == 0.0

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.repr()
    #-----------------------------------------------------------------
    assert elements.repr() == "1:H 2:C 3:O"

    #-----------------------------------------------------------------
    # Testing constructor ElementalComposition(list[tuple[Element, float]])
    #-----------------------------------------------------------------
    elements = ElementalComposition([
        ( Element("Ca"), 4.0 ),
        ( Element("Mg"), 5.0 ),
        ( Element("F"),  6.0 ),
    ])

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.symbols()
    #-----------------------------------------------------------------
    assert elements.symbols() == ["Ca", "Mg", "F"]

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.coefficients()
    #-----------------------------------------------------------------
    assert elements.coefficients() == [4.0, 5.0, 6.0]

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.coefficient()
    #-----------------------------------------------------------------
    assert elements.coefficient("Ca") == 4.0
    assert elements.coefficient("Mg") == 5.0
    assert elements.coefficient("F")  == 6.0
    assert elements.coefficient("X")  == 0.0

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.coefficient()
    #-----------------------------------------------------------------
    assert elements.coefficient("Ca") == 4.0
    assert elements.coefficient("Mg") == 5.0
    assert elements.coefficient("F")  == 6.0
    assert elements.coefficient("X")  == 0.0

    #-----------------------------------------------------------------
    # Testing method ElementalComposition.__iter__
    #-----------------------------------------------------------------
    for element, coefficient in elements:
        assert elements.coefficient(element.symbol()) == coefficient
