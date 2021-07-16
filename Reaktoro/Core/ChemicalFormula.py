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


def testChemicalFormula():

    formula = ChemicalFormula("H2O")
    assert formula.str()             == "H2O"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 2
    assert formula.coefficient("H")  == 2
    assert formula.coefficient("O")  == 1
    assert formula.equivalent("H2O(aq)")
    assert formula.equivalent("H2O(l)")
    assert formula.equivalent("HHO")
    assert formula.equivalent("HOH")
    assert formula.equivalent("OH2")

    formula = ChemicalFormula("CaCO3")
    assert formula.str()             == "CaCO3"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 3
    assert formula.coefficient("C")  == 1
    assert formula.coefficient("Ca") == 1
    assert formula.coefficient("O")  == 3
    assert formula.equivalent("CaCOOO")
    assert formula.equivalent("CaOOOC")
    assert formula.equivalent("Ca(CO3)")
    assert formula.equivalent("Ca(CO3)(aq)")
    assert formula.equivalent("Ca(CO3)(s)")

    formula = ChemicalFormula("HCO3-")
    assert formula.str()             == "HCO3-"
    assert formula.charge()          == -1
    assert len(formula.elements())   == 3
    assert formula.coefficient("C")  == 1
    assert formula.coefficient("H")  == 1
    assert formula.coefficient("O")  == 3
    assert formula.equivalent("HCO3-(aq)")
    assert formula.equivalent("HCO3[-](aq)")
    assert formula.equivalent("HCO3[-]")
    assert formula.equivalent("HCOOO-")

    formula = ChemicalFormula("H+")
    assert formula.str()             == "H+"
    assert formula.charge()          == 1
    assert len(formula.elements())   == 1
    assert formula.coefficient("H")  == 1
    assert formula.equivalent("H+(aq)")
    assert formula.equivalent("H[+]")

    formula = ChemicalFormula("e-")
    assert formula.str()             == "e-"
    assert formula.charge()          == -1
    assert len(formula.elements())   == 0
    assert formula.equivalent("e-(aq)")
    assert formula.equivalent("e[-]")

    formula = ChemicalFormula("Na+")
    assert formula.str()             == "Na+"
    assert formula.charge()          == 1
    assert len(formula.elements())   == 1
    assert formula.coefficient("Na") == 1
    assert formula.equivalent("Na+(aq)")
    assert formula.equivalent("Na+(pl)")
    assert formula.equivalent("Na[+](aq)")
    assert formula.equivalent("Na[+]")

    formula = ChemicalFormula("Cl-")
    assert formula.str()             == "Cl-"
    assert formula.charge()          == -1
    assert len(formula.elements())   == 1
    assert formula.coefficient("Cl") == 1
    assert formula.equivalent("Cl[-](aq)")
    assert formula.equivalent("Cl-(aq)")
    assert formula.equivalent("Cl-(pl)")

    formula = ChemicalFormula("CO3--")
    assert formula.str()             == "CO3--"
    assert formula.charge()          == -2
    assert len(formula.elements())   == 2
    assert formula.coefficient("C")  == 1
    assert formula.coefficient("O")  == 3
    assert formula.equivalent("CO3-2(aq)")
    assert formula.equivalent("CO3-2")
    assert formula.equivalent("CO3[2-](aq)")

    formula = ChemicalFormula("Fe+++")
    assert formula.str()             == "Fe+++"
    assert formula.charge()          == 3
    assert len(formula.elements())   == 1
    assert formula.coefficient("Fe") == 1
    assert formula.equivalent("Fe+++(aq)")
    assert formula.equivalent("Fe+3(aq)")
    assert formula.equivalent("Fe[3+](aq)")
    assert formula.equivalent("Fe+3")
    assert formula.equivalent("Fe[3+]")

    formula = ChemicalFormula("(CaMg)(CO3)2")
    assert formula.str()             == "(CaMg)(CO3)2"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 4
    assert formula.coefficient("C")  == 2
    assert formula.coefficient("Ca") == 1
    assert formula.coefficient("Mg") == 1
    assert formula.coefficient("O")  == 6
    assert formula.equivalent("CaMg(CO3)2(s)")
    assert formula.equivalent("CaMg(CO3)2")
    assert formula.equivalent("CaMgCO3CO3(s)")
    assert formula.equivalent("CaMgCO3CO3")

    formula = ChemicalFormula("CH3COOH")
    assert formula.str()             == "CH3COOH"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 3
    assert formula.coefficient("C")  == 2
    assert formula.coefficient("H")  == 4
    assert formula.coefficient("O")  == 2
    assert formula.equivalent("CH3COOH(aq)")
    assert formula.equivalent("CHHHCOOH(aq)")
    assert formula.equivalent("CHHHCOOH")

    formula = ChemicalFormula("Al2.5Si0.5O4.75")
    assert formula.str()             == "Al2.5Si0.5O4.75"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 3
    assert formula.coefficient("Al") == 2.5
    assert formula.coefficient("Si") == 0.5
    assert formula.coefficient("O")  == 4.75
    assert formula.equivalent("Al2.5Si0.5O4.75(s)")
    assert formula.equivalent("Si0.5Al2.5O4.75(s)")

    formula = ChemicalFormula("Fe4Al18Si7.5O48H4")
    assert formula.str()             == "Fe4Al18Si7.5O48H4"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 5
    assert formula.coefficient("Fe") == 4
    assert formula.coefficient("Al") == 18
    assert formula.coefficient("Si") == 7.5
    assert formula.coefficient("O")  == 48
    assert formula.coefficient("H")  == 4
    assert formula.equivalent("Fe4Al18Si7.5O48H4(s)")
    assert formula.equivalent("Al18Fe4H4Si7.5O48(s)")

    formula = ChemicalFormula("Mg4Al18Si7.5O48H4")
    assert formula.str()             == "Mg4Al18Si7.5O48H4"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 5
    assert formula.coefficient("Mg") == 4
    assert formula.coefficient("Al") == 18
    assert formula.coefficient("Si") == 7.5
    assert formula.coefficient("O")  == 48
    assert formula.coefficient("H")  == 4
    assert formula.equivalent("Mg4Al18Si7.5O48H4(s)")

    formula = ChemicalFormula("Mn4Al18Si7.5O48H4")
    assert formula.str()             == "Mn4Al18Si7.5O48H4"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 5
    assert formula.coefficient("Mn") == 4
    assert formula.coefficient("Al") == 18
    assert formula.coefficient("Si") == 7.5
    assert formula.coefficient("O")  == 48
    assert formula.coefficient("H")  == 4

    formula = ChemicalFormula("Ca0.5Al1Si2O6")
    assert formula.str()             == "Ca0.5Al1Si2O6"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 4
    assert formula.coefficient("Ca") == 0.5
    assert formula.coefficient("Al") == 1
    assert formula.coefficient("Si") == 2
    assert formula.coefficient("O")  == 6

    formula = ChemicalFormula("K0.5Fe5Al2Si8O30.5H12.5")
    assert formula.str()             == "K0.5Fe5Al2Si8O30.5H12.5"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 6
    assert formula.coefficient("K")  == 0.5
    assert formula.coefficient("Fe") == 5
    assert formula.coefficient("Al") == 2
    assert formula.coefficient("Si") == 8
    assert formula.coefficient("O")  == 30.5
    assert formula.coefficient("H")  == 12.5

    formula = ChemicalFormula("K0.5Mg5Al2Si8O30.5H12.5")
    assert formula.str()             == "K0.5Mg5Al2Si8O30.5H12.5"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 6
    assert formula.coefficient("K")  == 0.5
    assert formula.coefficient("Mg") == 5
    assert formula.coefficient("Al") == 2
    assert formula.coefficient("Si") == 8
    assert formula.coefficient("O")  == 30.5
    assert formula.coefficient("H")  == 12.5

    formula = ChemicalFormula("Mg3.5Al9Si1.5O20")
    assert formula.str()             == "Mg3.5Al9Si1.5O20"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 4
    assert formula.coefficient("Mg") == 3.5
    assert formula.coefficient("Al") == 9
    assert formula.coefficient("Si") == 1.5
    assert formula.coefficient("O")  == 20

    formula = ChemicalFormula("Fe3.5Al9Si1.5O20")
    assert formula.str()             == "Fe3.5Al9Si1.5O20"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 4
    assert formula.coefficient("Fe") == 3.5
    assert formula.coefficient("Al") == 9
    assert formula.coefficient("Si") == 1.5
    assert formula.coefficient("O")  == 20

    formula = ChemicalFormula("Fe0.875S1")
    assert formula.str()             == "Fe0.875S1"
    assert formula.charge()          == 0
    assert len(formula.elements())   == 2
    assert formula.coefficient("Fe") == 0.875
    assert formula.coefficient("S")  == 1

    assert ChemicalFormula.equivalent("Ca++", "Ca+2")
    assert ChemicalFormula.equivalent("Ca++", "Ca[2+]")
    assert ChemicalFormula.equivalent("Ca++", "Ca[2+](aq)")

    assert ChemicalFormula.equivalent("CO3--", "CO3-2")
    assert ChemicalFormula.equivalent("CO3--", "CO3[2-]")
    assert ChemicalFormula.equivalent("CO3--", "CO3[2-](aq)")

    assert ChemicalFormula.equivalent("Fe+++", "Fe+3")
    assert ChemicalFormula.equivalent("Fe+++", "Fe[3+]")
    assert ChemicalFormula.equivalent("Fe+++", "Fe[3+](aq)")

    assert ChemicalFormula.equivalent("H+", "H+1")
    assert ChemicalFormula.equivalent("H+", "H[+]")
    assert ChemicalFormula.equivalent("H+", "H[+](aq)")

    assert ChemicalFormula.equivalent("OH-", "OH-1")
    assert ChemicalFormula.equivalent("OH-", "OH[-]")
    assert ChemicalFormula.equivalent("OH-", "OH[-](aq)")

    assert ChemicalFormula.equivalent("CO2", "CO2(g)")
    assert ChemicalFormula.equivalent("CO2", "COO")
