#!/usr/bin/env python3

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


from parse_nasa_database import *


def test_formatNumber():
    assert formatNumber( 1.234) ==  1.234
    assert formatNumber( 1.000) ==  1
    assert formatNumber(-1.234) == -1.234
    assert formatNumber(-1.000) == -1
    assert type(formatNumber( 1.234)) == float
    assert type(formatNumber( 1.000)) == int
    assert type(formatNumber(-1.234)) == float
    assert type(formatNumber(-1.000)) == int


def test_correctNasaSpeciesName():
    assert correctNasaSpeciesName("NaCL") == "NaCl"
    assert correctNasaSpeciesName("WCL6") == "WCl6"
    assert correctNasaSpeciesName("AL(cr)") == "Al(cr)"
    assert correctNasaSpeciesName("ALBr3(cr)") == "AlBr3(cr)"


def test_identifyCommonSpeciesName():
    assert "BaCO3"     == identifyCommonSpeciesName(["BaCO3(a)", "BaCO3(b)", "BaCO3(c)", "BaCO3(L)"])
    assert "Cr2O3"     == identifyCommonSpeciesName(["Cr2O3(I')", "Cr2O3(I)", "Cr2O3(I)", "Cr2O3(I)", "Cr2O3(L)"])
    assert "Li3ALF6"   == identifyCommonSpeciesName(["Li3ALF6(IV)", "Li3ALF6(III)", "Li3ALF6(II)", "Li3ALF6(I)", "Li3ALF6(L)"])
    assert "Ag"        == identifyCommonSpeciesName(["Ag(cr)", "Ag(L)"])
    assert "Fe2O3(cr)" == identifyCommonSpeciesName(["Fe2O3(cr)", "Fe2O3(cr)"])
    assert "K2Si2O5"   == identifyCommonSpeciesName(["K2Si2O5(a)", "K2Si2O5(b)", "K2Si2O5(c)", "K2Si2O5(L)"])
    assert "SrCO3"     == identifyCommonSpeciesName(['SrCO3(a)', 'SrCO3(b)', 'SrCO3(c)', 'SrCO3(L)'])
    assert "SiO2"      == identifyCommonSpeciesName(['SiO2(a-qz)', 'SiO2(b-qz)', 'SiO2(b-crt)', 'SiO2(L)'])
    assert "K(HF2)"    == identifyCommonSpeciesName(['K(HF2)(a)', 'K(HF2)(b)', 'K(HF2)(L)'])


def test_determineSpeciesNameSuffix():
    pass


def test_combineSpeciesBlocksWhenPossible():
    pass


def test_segment():
    l = list(range(10))
    assert segment(l, 1, 4) == l[1:5]
    assert segment(l, 6, 2) == l[6:8]


def test_parseFortranScientificNumber():
    assert parseFortranScientificNumber("1.23")             ==  1.23
    assert parseFortranScientificNumber("1.23e5")           ==  1.23e5
    assert parseFortranScientificNumber("1.23e+4")          ==  1.23e+4
    assert parseFortranScientificNumber("1.23e-4")          ==  1.23e-4
    assert parseFortranScientificNumber("1.23014D+4")       ==  1.23014e+4
    assert parseFortranScientificNumber("1.23014D-4")       ==  1.23014e-4
    assert parseFortranScientificNumber(".23014D+4")        ==   .23014e+4
    assert parseFortranScientificNumber(".23014D-4")        ==   .23014e-4
    assert parseFortranScientificNumber("-3.309926370D+05") == -3.309926370e+05
    assert parseFortranScientificNumber(" 9.820086420D+02") ==  9.820086420e+02
    assert parseFortranScientificNumber(" 1.381179917D+00") ==  1.381179917e+00
    assert parseFortranScientificNumber(" 6.170899990D-04") ==  6.170899990e-04
    assert parseFortranScientificNumber("-1.688114600D-07") == -1.688114600e-07
    assert parseFortranScientificNumber("-3.309926370e+05") == -3.309926370e+05
    assert parseFortranScientificNumber(" 9.820086420e+02") ==  9.820086420e+02
    assert parseFortranScientificNumber(" 1.381179917e+00") ==  1.381179917e+00
    assert parseFortranScientificNumber(" 6.170899990e-04") ==  6.170899990e-04
    assert parseFortranScientificNumber("-1.688114600e-07") == -1.688114600e-07


numbers = "8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06"

def test_getStringBetweenColumns():
    assert getStringBetweenColumns(numbers,  1, 15) == "8.906039290D+04"
    assert getStringBetweenColumns(numbers, 16, 31) == "-9.750803930D+02"
    assert getStringBetweenColumns(numbers, 33, 33) == "5"


def test_getIntegerBetweenColumns():
    assert getIntBetweenColumns(numbers, 3, 4) == 90
    assert getIntBetweenColumns(numbers, 1, 1) == 8


def test_getDoubleBetweenColumns():
    assert getFloatBetweenColumns(numbers,  1, 15) ==  8.906039290e+04
    assert getFloatBetweenColumns(numbers, 16, 31) == -9.750803930e+02
    assert getFloatBetweenColumns(numbers, 33, 33) ==  5.0


def test_isCommentLine():
    assert isCommentLine("! Some text")
    assert isCommentLine("# Some text")
    assert isCommentLine("    ! Some text not trimmed ")
    assert isCommentLine("    # Some text not trimmed ")


def test_parseFormula():
    assert parseFormula("E   1.00    0.00    0.00    0.00    0.00") == [("E", 1.0)]
    assert parseFormula("AL  1.00H   2.00F   1.00    0.00    0.00") == [("AL", 1.0), ("H", 2.0), ("F", 1.0)]
    assert parseFormula("B   3.00N   3.00H   6.00    0.00    0.00") == [("B", 3.0), ("N", 3.0), ("H", 6.0)]
    assert parseFormula("NA  1.00O   1.00H   1.00E  -1.00    0.00") == [("NA", 1.0), ("O", 1.0), ("H", 1.0), ("E", -1.0)]
    assert parseFormula("ZR  1.00O   1.00E  -1.00    0.00    0.00") == [("ZR", 1.0), ("O", 1.0), ("E", -1.0)]
    assert parseFormula("Xx 12.00Ab  1.00Z   9.00Yy  1.20Qq  3.45") == [("Xx", 12.0), ("Ab", 1.0), ("Z", 9.0), ("Yy", 1.2), ("Qq", 3.45)]
    assert parseFormula("N 1.5617O .41959AR.00937C .00032  .00000") == [("N", 1.5617), ("O", 0.41959), ("AR", 0.00937), ("C", 0.00032)]
    assert parseFormula("C  12.00H  23.00    0.00    0.00    0.00") == [("C", 12.0), ("H", 23.0)]


def test_identifyAggregateState():
    assert identifyAggregateState("O2", 0) == "Gas"
    assert identifyAggregateState("H2", 0) == "Gas"

    assert identifyAggregateState("O2(L)", 0)  == "Gas"   # the aggregate code zero implies gas no matter how the species is name
    assert identifyAggregateState("Mg(cr)", 0) == "Gas"   # the aggregate code zero implies gas no matter how the species is name

    assert identifyAggregateState("Mg(L)", 1) == "Liquid"
    assert identifyAggregateState("Mg(l)", 1) == "Liquid"
    assert identifyAggregateState("C2H2(L),acetyle", 4) == "Liquid"

    assert identifyAggregateState("Mg(cr)", 1) == "Solid"
    assert identifyAggregateState("MgO(s)", 1) == "Solid"
    assert identifyAggregateState("CaSO4(II)", 3) == "Solid"
    assert identifyAggregateState("CaSO4(IV)", 3) == "Solid"
    assert identifyAggregateState("CsCL(a)", 1) == "Solid"
    assert identifyAggregateState("CsCL(b)", 2) == "Solid"

    assert identifyAggregateState("RP-1", 2) == "Liquid"
    assert identifyAggregateState("JP-4", 3) == "Liquid"
    assert identifyAggregateState("JP-5", 4) == "Liquid"
    assert identifyAggregateState("IRFNA", 5) == "Liquid"


def test_createNasaPolynomial():

    lines00 = [
        "    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10764.801",
        " 8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06",
        " 7.219318660D-10 7.329846740D-14                -3.425302600D+04-9.387539141D+00"
    ]

    params = createNasaPolynomial(lines00)

    assert params.Tmin ==  200.0
    assert params.Tmax ==  1000.0
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  8.906039290e+04
    assert params.a2   == -9.750803930e+02
    assert params.a3   ==  5.870772320e+00
    assert params.a4   ==  6.896663690e-03
    assert params.a5   == -4.045817280e-06
    assert params.a6   ==  7.219318660e-10
    assert params.a7   ==  7.329846740e-14
    assert params.b1   == -3.425302600e+04
    assert params.b2   == -9.387539141e+00

    lines01 = [
        "    287.700    288.5007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-1.494237648D+11 0.000000000D+00 8.652205900D+06-3.521604880D+04 3.968501633D+01",
        " 0.000000000D+00 0.000000000D+00                -1.866195945D+09-4.140143280D+07"
    ]

    params = createNasaPolynomial(lines01)

    assert params.Tmin ==  287.700
    assert params.Tmax ==  288.500
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   == -1.494237648e+11
    assert params.a2   ==  0.000000000e+00
    assert params.a3   ==  8.652205900e+06
    assert params.a4   == -3.521604880e+04
    assert params.a5   ==  3.968501633e+01
    assert params.a6   ==  0.000000000e+00
    assert params.a7   ==  0.000000000e+00
    assert params.b1   == -1.866195945e+09
    assert params.b2   == -4.140143280e+07

    lines02 = [
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10772.157",
        "-3.834452340D+04-9.563261140D+02 8.200962780D+00-2.062130802D-04 9.872288990D-09",
        " 8.158366760D-12-7.527519660D-16                -3.423564020D+04-2.224772278D+01"
    ]

    params = createNasaPolynomial(lines02)

    assert params.Tmin ==  1000.000
    assert params.Tmax ==  6000.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   == -3.834452340e+04
    assert params.a2   == -9.563261140e+02
    assert params.a3   ==  8.200962780e+00
    assert params.a4   == -2.062130802e-04
    assert params.a5   ==  9.872288990e-09
    assert params.a6   ==  8.158366760e-12
    assert params.a7   == -7.527519660e-16
    assert params.b1   == -3.423564020e+04
    assert params.b2   == -2.224772278e+01

    lines03 = [
        "   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428",
        " 2.291781634D+09-1.608862960D+06 4.312466360D+02-5.396508990D-02 3.531856210D-06",
        "-1.164403850D-10 1.527134223D-15                 1.258651434D+07-3.692101610D+03"
    ]

    params = createNasaPolynomial(lines03)

    assert params.Tmin ==  6000.000
    assert params.Tmax ==  20000.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  2.291781634e+09
    assert params.a2   == -1.608862960e+06
    assert params.a3   ==  4.312466360e+02
    assert params.a4   == -5.396508990e-02
    assert params.a5   ==  3.531856210e-06
    assert params.a6   == -1.164403850e-10
    assert params.a7   ==  1.527134223e-15
    assert params.b1   ==  1.258651434e+07
    assert params.b2   == -3.692101610e+03

    lines04 = [
        "    298.150    700.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 5.676714650D+06-6.561919460D+04 2.971622148D+02-6.132148300D-01 6.373067250D-04",
        "-2.569754132D-07 0.000000000D+00                 3.113488256D+05-1.716694497D+03"
    ]

    params = createNasaPolynomial(lines04)

    assert params.Tmin ==  298.150
    assert params.Tmax ==  700.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  5.676714650e+06
    assert params.a2   == -6.561919460e+04
    assert params.a3   ==  2.971622148e+02
    assert params.a4   == -6.132148300e-01
    assert params.a5   ==  6.373067250e-04
    assert params.a6   == -2.569754132e-07
    assert params.a7   ==  0.000000000e+00
    assert params.b1   ==  3.113488256e+05
    assert params.b2   == -1.716694497e+03

    lines05 = [
        "    700.000   1500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 8.246680410D+05-1.965549089D+03 5.659967790D+00 1.366241998D-02-1.241835551D-05",
        " 6.313550760D-09-1.192405400D-12                -3.594772610D+03-3.071970867D+01"
    ]

    params = createNasaPolynomial(lines05)

    assert params.Tmin ==  700.000
    assert params.Tmax ==  1500.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  8.246680410e+05
    assert params.a2   == -1.965549089e+03
    assert params.a3   ==  5.659967790e+00
    assert params.a4   ==  1.366241998e-02
    assert params.a5   == -1.241835551e-05
    assert params.a6   ==  6.313550760e-09
    assert params.a7   == -1.192405400e-12
    assert params.b1   == -3.594772610e+03
    assert params.b2   == -3.071970867e+01

    lines06 = [
        "   1500.000   2500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 4.325258320D+09-1.355254539D+07 1.757503878D+04-1.205560629D+01 4.623646210D-03",
        "-9.389890050D-07 7.891119490D-11                 8.506076540D+07-1.227417790D+05"
    ]

    params = createNasaPolynomial(lines06)

    assert params.Tmin ==  1500.000
    assert params.Tmax ==  2500.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  4.325258320e+09
    assert params.a2   == -1.355254539e+07
    assert params.a3   ==  1.757503878e+04
    assert params.a4   == -1.205560629e+01
    assert params.a5   ==  4.623646210e-03
    assert params.a6   == -9.389890050e-07
    assert params.a7   ==  7.891119490e-11
    assert params.b1   ==  8.506076540e+07
    assert params.b2   == -1.227417790e+05

    lines07 = [
        "    200.000    306.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-2.112979833D+06 0.000000000D+00 2.207431906D+02-1.164973854D+00 1.853971249D-03",
        " 0.000000000D+00 0.000000000D+00                -1.746831159D+05-9.949025210D+02"
    ]

    params = createNasaPolynomial(lines07)

    assert params.Tmin ==  200.000
    assert params.Tmax ==  306.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   == -2.112979833e+06
    assert params.a2   ==  0.000000000e+00
    assert params.a3   ==  2.207431906e+02
    assert params.a4   == -1.164973854e+00
    assert params.a5   ==  1.853971249e-03
    assert params.a6   ==  0.000000000e+00
    assert params.a7   ==  0.000000000e+00
    assert params.b1   == -1.746831159e+05
    assert params.b2   == -9.949025210e+02

    lines08 = [
        "    306.000    310.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 6.705915562D+03-4.303760534D+01 6.919229155D-02",
        " 0.000000000D+00 0.000000000D+00                -8.349875160D+05-2.844167579D+04"
    ]

    params = createNasaPolynomial(lines08)

    assert params.Tmin ==  306.000
    assert params.Tmax ==  310.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  0.000000000e+00
    assert params.a2   ==  0.000000000e+00
    assert params.a3   ==  6.705915562e+03
    assert params.a4   == -4.303760534e+01
    assert params.a5   ==  6.919229155e-02
    assert params.a6   ==  0.000000000e+00
    assert params.a7   ==  0.000000000e+00
    assert params.b1   == -8.349875160e+05
    assert params.b2   == -2.844167579e+04

    lines09 = [
        "    310.000    335.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 2.443570337D+02-1.399445548D+00 2.113509996D-03",
        " 0.000000000D+00 0.000000000D+00                -1.665032895D+05-1.059172219D+03"
    ]

    params = createNasaPolynomial(lines09)

    assert params.Tmin ==  310.000
    assert params.Tmax ==  335.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   ==  0.000000000e+00
    assert params.a2   ==  0.000000000e+00
    assert params.a3   ==  2.443570337e+02
    assert params.a4   == -1.399445548e+00
    assert params.a5   ==  2.113509996e-03
    assert params.a6   ==  0.000000000e+00
    assert params.a7   ==  0.000000000e+00
    assert params.b1   == -1.665032895e+05
    assert params.b2   == -1.059172219e+03

    lines10 = [
        "    335.000   2705.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-3.415474875D+05 0.000000000D+00 1.616932327D+01-1.517828471D-03 1.014852348D-06",
        " 0.000000000D+00 0.000000000D+00                -1.430478214D+05-8.374919535D+01"
    ]

    params = createNasaPolynomial(lines10)

    assert params.Tmin ==  335.000
    assert params.Tmax ==  2705.000
    assert params.qN   ==  7
    assert params.q1   == -2.0
    assert params.q2   == -1.0
    assert params.q3   ==  0.0
    assert params.q4   ==  1.0
    assert params.q5   ==  2.0
    assert params.q6   ==  3.0
    assert params.q7   ==  4.0
    assert params.a1   == -3.415474875e+05
    assert params.a2   ==  0.000000000e+00
    assert params.a3   ==  1.616932327e+01
    assert params.a4   == -1.517828471e-03
    assert params.a5   ==  1.014852348e-06
    assert params.a6   ==  0.000000000e+00
    assert params.a7   ==  0.000000000e+00
    assert params.b1   == -1.430478214e+05
    assert params.b2   == -8.374919535e+01


def test_createNasaSpecies():
    lines11 = [
        "MgO               Gurvich,1996a pt1 p398 pt2 p318.",
        " 3 tpis96 MG  1.00O   1.00    0.00    0.00    0.00 0   40.3044000      32261.307",
        "    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8909.107",
        " 3.513659740D+05-5.287197160D+03 3.382060060D+01-8.400489630D-02 1.210016160D-04",
        "-7.630795020D-08 1.701022862D-11                 2.790679519D+04-1.624886199D+02",
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8909.107",
        "-1.586738367D+07 3.420468100D+04-1.774087677D+01 7.004963050D-03-1.104138249D-06",
        " 8.957488530D-11-3.052513649D-15                -2.300504434D+05 1.738984472D+02",
        "   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8909.107",
        " 2.290059050D+06-2.073499632D+04 1.444150005D+01-1.490609900D-03 1.052119343D-07",
        "-3.523030610D-12 4.613111760D-17                 1.490218815D+05-8.007281730D+01",
    ]

    species = createNasaSpecies(lines11)

    assert species.name == "MgO"
    assert species.comment == "Gurvich,1996a pt1 p398 pt2 p318."
    assert species.idcode == "tpis96"
    assert species.elements == [("MG", 1.0), ("O", 1.0)]
    assert species.aggregatecode == 0
    assert species.aggregatestate == "Gas"
    assert species.speciestype == ""
    assert species.molarmass == 40.3044000
    assert species.dHf == 32261.307
    assert species.dH0 == 8909.107
    assert species.H0 == 0.0
    assert len(species.polynomials) == 3
    assert repr(species.polynomials[0]) == repr(createNasaPolynomial(segment(lines11, 2, 3)))
    assert repr(species.polynomials[1]) == repr(createNasaPolynomial(segment(lines11, 5, 3)))
    assert repr(species.polynomials[2]) == repr(createNasaPolynomial(segment(lines11, 8, 3)))

    lines12 = [
        "Mg(cr)            Hexagonal. Ref-Elm. Alcock,1993.",
        " 2 srd 93 MG  1.00    0.00    0.00    0.00    0.00 1   24.3050000          0.000",
        "    100.000    298.1507 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4979.161",
        "-5.412225134D+03 0.000000000D+00 1.458173723D+00 1.330204666D-02-4.098858502D-05",
        " 4.754339101D-08 0.000000000D+00                -7.759472010D+02-6.989702348D+00",
        "    298.150    923.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4979.161",
        "-2.860060304D+04 0.000000000D+00 3.398877384D+00-7.243962663D-04 1.405254188D-06",
        " 0.000000000D+00 0.000000000D+00                -1.089519906D+03-1.545973664D+01",
    ]

    species = createNasaSpecies(lines12)

    assert species.name == "Mg(cr)"
    assert species.comment == "Hexagonal. Ref-Elm. Alcock,1993."
    assert species.idcode == "srd 93"
    assert species.elements == [("MG", 1.0)]
    assert species.aggregatecode == 1
    assert species.aggregatestate == "Solid"
    assert species.speciestype == ""
    assert species.molarmass == 24.3050000
    assert species.dHf == 0.000
    assert species.dH0 == 4979.161
    assert species.H0 == 0.0
    assert len(species.polynomials) == 2
    assert repr(species.polynomials[0]) == repr(createNasaPolynomial(segment(lines12, 2, 3)))
    assert repr(species.polynomials[1]) == repr(createNasaPolynomial(segment(lines12, 5, 3)))

    lines13 = [
        "Jet-A(g)          McBride,1996. Faith,1971. Gracia-Salcedo,1988.          React.",
        " 2 g 8/01 C  12.00H  23.00    0.00    0.00    0.00 0  167.3110200    -249657.000",
        "    273.150   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        "-6.068695590D+05 8.328259590D+03-4.312321270D+01 2.572390455D-01-2.629316040D-04",
        " 1.644988940D-07-4.645335140D-11                -7.606962760D+04 2.794305937D+02",
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 1.858356102D+07-7.677219890D+04 1.419826133D+02-7.437524530D-03 5.856202550D-07",
        " 1.223955647D-11-3.149201922D-15                 4.221989520D+05-8.986061040D+02",
    ]

    species = createNasaSpecies(lines13)

    assert species.name == "Jet-A(g)"
    assert species.comment == "McBride,1996. Faith,1971. Gracia-Salcedo,1988.          React."
    assert species.idcode == "g 8/01"
    assert species.elements == [("C", 12.0), ("H", 23.0)]
    assert species.aggregatecode == 0
    assert species.aggregatestate == "Gas"
    assert species.speciestype == ""
    assert species.molarmass == 167.3110200
    assert species.dHf == -249657.000
    assert species.dH0 == 0.0
    assert species.H0 == 0.0
    assert len(species.polynomials) == 2
    assert repr(species.polynomials[0]) == repr(createNasaPolynomial(segment(lines13, 2, 3)))
    assert repr(species.polynomials[1]) == repr(createNasaPolynomial(segment(lines13, 5, 3)))

    lines14 = [
        "NaCN(II)          Lambda trans@288.5K. Chase,1998(3/66) pp631-3. Messer,1941.",
        " 6 g 8/01 NA  1.00C   1.00N   1.00    0.00    0.00 1   49.0071700     -90709.000",
        "    197.700    245.9007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        " 4.995073610D+08-1.007656780D+07 7.458386410D+04-2.134406974D+02-6.965856220D-02",
        " 1.540182072D-03-2.204709373D-06                 4.493712560D+07-3.984423000D+05",
        "    245.900    273.1007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-4.546907220D+08 4.339535660D+06 5.476070770D+02-1.324240852D+02 5.045446840D-01",
        "-5.791737950D-04 0.000000000D+00                -2.385188911D+07 3.106410058D+04",
        "    273.100    284.2007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-4.415096320D+08 0.000000000D+00 3.452631670D+04-1.661761169D+02 2.250608921D-01",
        " 0.000000000D+00 0.000000000D+00                -6.388306670D+06-1.596447272D+05",
        "    284.200    286.3007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        " 1.647379186D+09 0.000000000D+00-6.100370390D+04 1.429224260D+02 0.000000000D+00",
        " 0.000000000D+00 0.000000000D+00                 1.735057318D+07 3.142435378D+05",
        "    286.300    287.7007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        " 1.607424258D+11 0.000000000D+00-8.629992400D+06 3.290155880D+04-3.355904228D+01",
        " 0.000000000D+00 0.000000000D+00                 1.946284579D+09 4.175641230D+07",
        "    287.700    288.5007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-1.494237648D+11 0.000000000D+00 8.652205900D+06-3.521604880D+04 3.968501633D+01",
        " 0.000000000D+00 0.000000000D+00                -1.866195945D+09-4.140143280D+07",
    ]

    species = createNasaSpecies(lines14)

    assert species.name == "NaCN(II)"
    assert species.comment == "Lambda trans@288.5K. Chase,1998(3/66) pp631-3. Messer,1941."
    assert species.idcode == "g 8/01"
    assert species.elements == [("NA", 1.0), ("C", 1.0), ("N", 1.0)]
    assert species.aggregatecode == 1
    assert species.aggregatestate == "Solid"
    assert species.speciestype == ""
    assert species.molarmass == 49.0071700
    assert species.dHf == -90709.000
    assert species.dH0 == 19422.128
    assert species.H0 == 0.0
    assert len(species.polynomials) == 6
    assert repr(species.polynomials[0]) == repr(createNasaPolynomial(segment(lines14, 2, 3)))
    assert repr(species.polynomials[1]) == repr(createNasaPolynomial(segment(lines14, 5, 3)))
    assert repr(species.polynomials[2]) == repr(createNasaPolynomial(segment(lines14, 8, 3)))
    assert repr(species.polynomials[3]) == repr(createNasaPolynomial(segment(lines14, 11, 3)))
    assert repr(species.polynomials[4]) == repr(createNasaPolynomial(segment(lines14, 14, 3)))
    assert repr(species.polynomials[5]) == repr(createNasaPolynomial(segment(lines14, 17, 3)))

    lines15 = [
        "C2H2(L),acetyle   Acetylene. McBride,1996 pp84,92.",
        " 0 g 6/96 C   2.00H   2.00    0.00    0.00    0.00 1   26.0372800     207599.000",
        "    192.350      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
    ]

    species = createNasaSpecies(lines15)

    assert species.name == "C2H2(l),acetyle"
    assert species.comment == "Acetylene. McBride,1996 pp84,92."
    assert species.idcode == "g 6/96"
    assert species.elements == [("C", 2.0), ("H", 2.0)]
    assert species.aggregatecode == 1
    assert species.aggregatestate == "Liquid"
    assert species.speciestype == ""
    assert species.molarmass == 26.0372800
    assert species.dHf == 0.0
    assert species.dH0 == 0.0
    assert species.H0 == 207599.000
    assert species.polynomials == []

    lines16 = [
        "N2O4(L)           Dinitrogen tetroxide. McBride,1996 pp85,93.",
        " 0 g 6/96 N   2.00O   4.00    0.00    0.00    0.00 1   92.0110000     -17549.000",
        "    298.150      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
    ]

    species = createNasaSpecies(lines16)

    assert species.name == "N2O4(l)"
    assert species.comment == "Dinitrogen tetroxide. McBride,1996 pp85,93."
    assert species.idcode == "g 6/96"
    assert species.elements == [("N", 2.0), ("O", 4.0)]
    assert species.aggregatecode == 1
    assert species.aggregatestate == "Liquid"
    assert species.speciestype == ""
    assert species.molarmass == 92.0110000
    assert species.dHf == 0.0
    assert species.dH0 == 0.0
    assert species.H0 == -17549.000
    assert species.polynomials == []


def test_createNasaSpeciesList():
    pass


def test_createFormulaString():
    assert "e-"      == createFormulaString([("E", 1)])
    assert "H2O"     == createFormulaString([("H", 2), ("O", 1)])
    assert "H+"      == createFormulaString([("H", 1), ("E", -1)])
    assert "CO3-2"   == createFormulaString([("C", 1), ("O", 3), ("E", 2)])
    assert "Mg"      == createFormulaString([("Mg", 1)])
    assert "Fe+3"    == createFormulaString([("Fe", 1), ("E", -3)])


def test_createElementsString():
    assert ""        == createElementsString([("E", 1)])
    assert "2:H 1:O" == createElementsString([("H", 2), ("O", 1)])
    assert "1:H"     == createElementsString([("H", 1), ("E", -1)])
    assert "1:C 3:O" == createElementsString([("C", 1), ("O", 3), ("E", 2)])
    assert "1:Mg"    == createElementsString([("Mg", 1)])
    assert "1:Fe"    == createElementsString([("Fe", 1), ("E", -3)])


def test_getCharge():
    assert -1 == getCharge([("E", 1)])
    assert  0 == getCharge([("H", 2), ("O", 1)])
    assert  1 == getCharge([("H", 1), ("E", -1)])
    assert -2 == getCharge([("C", 1), ("O", 3), ("E", 2)])
    assert  0 == getCharge([("Mg", 1)])
    assert  3 == getCharge([("Fe", 1), ("E", -3)])


def test_getTags():
    assert getTags("") == "product reactant"
    assert getTags("product") == "product"
    assert getTags("reactant") == "reactant"


def test_getStandardThermoModel():
    species = NasaSpecies()
    species.name = "Fe2O3(cr)"
    species.formula = "Fe2O3"
    species.elements = [("Fe", 2), ("O", 3)]
    species.speciestype = "product"
    species.aggregatestate = "Solid"
    species.dHf = -824248.0
    species.dH0 = 0.0

    poly0 = NasaPolynomial()
    poly0.Tmin = 960.0
    poly0.Tmax = 1800.0
    poly0.a1 = -28440090.66
    poly0.a2 = 125378.645
    poly0.a3 = -211.9080355
    poly0.a4 = 0.2221481558
    poly0.a5 = -0.0001206765164
    poly0.a6 = 3.48073805e-08
    poly0.a7 = -4.16529786e-12
    poly0.b1 = -848392.652
    poly0.b2 = 1433.002161

    poly1 = NasaPolynomial()
    poly1.Tmin = 1800.0
    poly1.Tmax = 6000.0
    poly1.a1 = 0.0
    poly1.a2 = 0.0
    poly1.a3 = 17.1143988
    poly1.a4 = 0.0
    poly1.a5 = 0.0
    poly1.a6 = 0.0
    poly1.a7 = 0.0
    poly1.b1 = -104159.9645
    poly1.b2 = -87.80613745

    species.polynomials = [poly0, poly1]

    data = getStandardThermoModel([species])

    assert data["Nasa"]["Polynomials"][0]["Tmin"] == poly0.Tmin
    assert data["Nasa"]["Polynomials"][0]["Tmax"] == poly0.Tmax
    assert data["Nasa"]["Polynomials"][0]["a1"] == poly0.a1
    assert data["Nasa"]["Polynomials"][0]["a2"] == poly0.a2
    assert data["Nasa"]["Polynomials"][0]["a3"] == poly0.a3
    assert data["Nasa"]["Polynomials"][0]["a4"] == poly0.a4
    assert data["Nasa"]["Polynomials"][0]["a5"] == poly0.a5
    assert data["Nasa"]["Polynomials"][0]["a6"] == poly0.a6
    assert data["Nasa"]["Polynomials"][0]["a7"] == poly0.a7
    assert data["Nasa"]["Polynomials"][0]["b1"] == poly0.b1
    assert data["Nasa"]["Polynomials"][0]["b2"] == poly0.b2

    assert data["Nasa"]["Polynomials"][1]["Tmin"] == poly1.Tmin
    assert data["Nasa"]["Polynomials"][1]["Tmax"] == poly1.Tmax
    assert data["Nasa"]["Polynomials"][1]["a1"] == poly1.a1
    assert data["Nasa"]["Polynomials"][1]["a2"] == poly1.a2
    assert data["Nasa"]["Polynomials"][1]["a3"] == poly1.a3
    assert data["Nasa"]["Polynomials"][1]["a4"] == poly1.a4
    assert data["Nasa"]["Polynomials"][1]["a5"] == poly1.a5
    assert data["Nasa"]["Polynomials"][1]["a6"] == poly1.a6
    assert data["Nasa"]["Polynomials"][1]["a7"] == poly1.a7
    assert data["Nasa"]["Polynomials"][1]["b1"] == poly1.b1
    assert data["Nasa"]["Polynomials"][1]["b2"] == poly1.b2
