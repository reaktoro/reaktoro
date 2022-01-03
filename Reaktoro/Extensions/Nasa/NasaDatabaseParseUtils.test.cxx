// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Extensions/Nasa/NasaDatabaseParseUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing NasaDatabaseParseUtils module", "[NasaDatabaseParseUtils]")
{
    using namespace NasaUtils;

    CHECK( convertFortranScientificNumberToDouble("-3.309926370D+05") == Approx(-3.309926370e+05) );
    CHECK( convertFortranScientificNumberToDouble(" 9.820086420D+02") == Approx( 9.820086420e+02) );
    CHECK( convertFortranScientificNumberToDouble(" 1.381179917D+00") == Approx( 1.381179917e+00) );
    CHECK( convertFortranScientificNumberToDouble(" 6.170899990D-04") == Approx( 6.170899990e-04) );
    CHECK( convertFortranScientificNumberToDouble("-1.688114600D-07") == Approx(-1.688114600e-07) );
    CHECK( convertFortranScientificNumberToDouble("-3.309926370e+05") == Approx(-3.309926370e+05) );
    CHECK( convertFortranScientificNumberToDouble(" 9.820086420e+02") == Approx( 9.820086420e+02) );
    CHECK( convertFortranScientificNumberToDouble(" 1.381179917e+00") == Approx( 1.381179917e+00) );
    CHECK( convertFortranScientificNumberToDouble(" 6.170899990e-04") == Approx( 6.170899990e-04) );
    CHECK( convertFortranScientificNumberToDouble("-1.688114600e-07") == Approx(-1.688114600e-07) );

    const auto numbers = "8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06";

    CHECK( getStringBetweenColumns(numbers,  1, 15) == "8.906039290D+04"  );
    CHECK( getStringBetweenColumns(numbers, 16, 31) == "-9.750803930D+02" );
    CHECK( getStringBetweenColumns(numbers, 33, 33) == "5"                );

    CHECK( getIntegerBetweenColumns(numbers, 3, 4) == 90 );
    CHECK( getIntegerBetweenColumns(numbers, 1, 1) == 8  );

    CHECK( getDoubleBetweenColumns(numbers,  1, 15) == Approx(8.906039290e+04)  );
    CHECK( getDoubleBetweenColumns(numbers, 16, 31) == Approx(-9.750803930e+02) );
    CHECK( getDoubleBetweenColumns(numbers, 33, 33) == Approx(5.0)              );

    CHECK( isCommentLine("! Some text") );
    CHECK( isCommentLine("# Some text") );
    CHECK( isCommentLine("    ! Some text not trimmed ") );
    CHECK( isCommentLine("    # Some text not trimmed ") );

    CHECK( parseFormula("E   1.00    0.00    0.00    0.00    0.00")
        == Pairs<String, double>{{"E", 1.0}} );

    CHECK( parseFormula("AL  1.00H   2.00F   1.00    0.00    0.00")
        == Pairs<String, double>{{"AL", 1.0}, {"H", 2.0}, {"F", 1.0}} );

    CHECK( parseFormula("B   3.00N   3.00H   6.00    0.00    0.00")
        == Pairs<String, double>{{"B", 3.0}, {"N", 3.0}, {"H", 6.0}} );

    CHECK( parseFormula("NA  1.00O   1.00H   1.00E  -1.00    0.00")
        == Pairs<String, double>{{"NA", 1.0}, {"O", 1.0}, {"H", 1.0}, {"E", -1.0}} );

    CHECK( parseFormula("ZR  1.00O   1.00E  -1.00    0.00    0.00")
        == Pairs<String, double>{{"ZR", 1.0}, {"O", 1.0}, {"E", -1.0}} );

    CHECK( parseFormula("Xx 12.00Ab  1.00Z   9.00Yy  1.20Qq  3.45")
        == Pairs<String, double>{{"Xx", 12.0}, {"Ab", 1.0}, {"Z", 9.0}, {"Yy", 1.2}, {"Qq", 3.45}} );

    CHECK( parseFormula("N 1.5617O .41959AR.00937C .00032  .00000")
        == Pairs<String, double>{{"N", 1.5617}, {"O", 0.41959}, {"AR", 0.00937}, {"C", 0.00032}} );

    CHECK( parseFormula("C  12.00H  23.00    0.00    0.00    0.00")
        == Pairs<String, double>{{"C", 12.0}, {"H", 23.0}} );

    Strings lines;
    NasaSpeciesThermoParams params;

    lines = {
        "    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10764.801",
        " 8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06",
        " 7.219318660D-10 7.329846740D-14                -3.425302600D+04-9.387539141D+00" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  200.0 );
    CHECK( params.Tmax ==  1000.0 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  8.906039290e+04 );
    CHECK( params.a2   == -9.750803930e+02 );
    CHECK( params.a3   ==  5.870772320e+00 );
    CHECK( params.a4   ==  6.896663690e-03 );
    CHECK( params.a5   == -4.045817280e-06 );
    CHECK( params.a6   ==  7.219318660e-10 );
    CHECK( params.a7   ==  7.329846740e-14 );
    CHECK( params.b1   == -3.425302600e+04 );
    CHECK( params.b2   == -9.387539141e+00 );

    lines = {
        "    287.700    288.5007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-1.494237648D+11 0.000000000D+00 8.652205900D+06-3.521604880D+04 3.968501633D+01",
        " 0.000000000D+00 0.000000000D+00                -1.866195945D+09-4.140143280D+07" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  287.700 );
    CHECK( params.Tmax ==  288.500 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   == -1.494237648e+11 );
    CHECK( params.a2   ==  0.000000000e+00 );
    CHECK( params.a3   ==  8.652205900e+06 );
    CHECK( params.a4   == -3.521604880e+04 );
    CHECK( params.a5   ==  3.968501633e+01 );
    CHECK( params.a6   ==  0.000000000e+00 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   == -1.866195945e+09 );
    CHECK( params.b2   == -4.140143280e+07 );

    lines = {
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10772.157",
        "-3.834452340D+04-9.563261140D+02 8.200962780D+00-2.062130802D-04 9.872288990D-09",
        " 8.158366760D-12-7.527519660D-16                -3.423564020D+04-2.224772278D+01" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  1000.000 );
    CHECK( params.Tmax ==  6000.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   == -3.834452340e+04 );
    CHECK( params.a2   == -9.563261140e+02 );
    CHECK( params.a3   ==  8.200962780e+00 );
    CHECK( params.a4   == -2.062130802e-04 );
    CHECK( params.a5   ==  9.872288990e-09 );
    CHECK( params.a6   ==  8.158366760e-12 );
    CHECK( params.a7   == -7.527519660e-16 );
    CHECK( params.b1   == -3.423564020e+04 );
    CHECK( params.b2   == -2.224772278e+01 );

    lines = {
        "   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428",
        " 2.291781634D+09-1.608862960D+06 4.312466360D+02-5.396508990D-02 3.531856210D-06",
        "-1.164403850D-10 1.527134223D-15                 1.258651434D+07-3.692101610D+03" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  6000.000 );
    CHECK( params.Tmax ==  20000.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  2.291781634e+09 );
    CHECK( params.a2   == -1.608862960e+06 );
    CHECK( params.a3   ==  4.312466360e+02 );
    CHECK( params.a4   == -5.396508990e-02 );
    CHECK( params.a5   ==  3.531856210e-06 );
    CHECK( params.a6   == -1.164403850e-10 );
    CHECK( params.a7   ==  1.527134223e-15 );
    CHECK( params.b1   ==  1.258651434e+07 );
    CHECK( params.b2   == -3.692101610e+03 );

    lines = {
        "    298.150    700.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 5.676714650D+06-6.561919460D+04 2.971622148D+02-6.132148300D-01 6.373067250D-04",
        "-2.569754132D-07 0.000000000D+00                 3.113488256D+05-1.716694497D+03" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  298.150 );
    CHECK( params.Tmax ==  700.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  5.676714650e+06 );
    CHECK( params.a2   == -6.561919460e+04 );
    CHECK( params.a3   ==  2.971622148e+02 );
    CHECK( params.a4   == -6.132148300e-01 );
    CHECK( params.a5   ==  6.373067250e-04 );
    CHECK( params.a6   == -2.569754132e-07 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   ==  3.113488256e+05 );
    CHECK( params.b2   == -1.716694497e+03 );

    lines = {
        "    700.000   1500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 8.246680410D+05-1.965549089D+03 5.659967790D+00 1.366241998D-02-1.241835551D-05",
        " 6.313550760D-09-1.192405400D-12                -3.594772610D+03-3.071970867D+01" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  700.000 );
    CHECK( params.Tmax ==  1500.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  8.246680410e+05 );
    CHECK( params.a2   == -1.965549089e+03 );
    CHECK( params.a3   ==  5.659967790e+00 );
    CHECK( params.a4   ==  1.366241998e-02 );
    CHECK( params.a5   == -1.241835551e-05 );
    CHECK( params.a6   ==  6.313550760e-09 );
    CHECK( params.a7   == -1.192405400e-12 );
    CHECK( params.b1   == -3.594772610e+03 );
    CHECK( params.b2   == -3.071970867e+01 );

    lines = {
        "   1500.000   2500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 4.325258320D+09-1.355254539D+07 1.757503878D+04-1.205560629D+01 4.623646210D-03",
        "-9.389890050D-07 7.891119490D-11                 8.506076540D+07-1.227417790D+05" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  1500.000 );
    CHECK( params.Tmax ==  2500.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  4.325258320e+09 );
    CHECK( params.a2   == -1.355254539e+07 );
    CHECK( params.a3   ==  1.757503878e+04 );
    CHECK( params.a4   == -1.205560629e+01 );
    CHECK( params.a5   ==  4.623646210e-03 );
    CHECK( params.a6   == -9.389890050e-07 );
    CHECK( params.a7   ==  7.891119490e-11 );
    CHECK( params.b1   ==  8.506076540e+07 );
    CHECK( params.b2   == -1.227417790e+05 );

    lines = {
        "    200.000    306.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-2.112979833D+06 0.000000000D+00 2.207431906D+02-1.164973854D+00 1.853971249D-03",
        " 0.000000000D+00 0.000000000D+00                -1.746831159D+05-9.949025210D+02" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  200.000 );
    CHECK( params.Tmax ==  306.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   == -2.112979833e+06 );
    CHECK( params.a2   ==  0.000000000e+00 );
    CHECK( params.a3   ==  2.207431906e+02 );
    CHECK( params.a4   == -1.164973854e+00 );
    CHECK( params.a5   ==  1.853971249e-03 );
    CHECK( params.a6   ==  0.000000000e+00 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   == -1.746831159e+05 );
    CHECK( params.b2   == -9.949025210e+02 );

    lines = {
        "    306.000    310.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 6.705915562D+03-4.303760534D+01 6.919229155D-02",
        " 0.000000000D+00 0.000000000D+00                -8.349875160D+05-2.844167579D+04" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  306.000 );
    CHECK( params.Tmax ==  310.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  0.000000000e+00 );
    CHECK( params.a2   ==  0.000000000e+00 );
    CHECK( params.a3   ==  6.705915562e+03 );
    CHECK( params.a4   == -4.303760534e+01 );
    CHECK( params.a5   ==  6.919229155e-02 );
    CHECK( params.a6   ==  0.000000000e+00 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   == -8.349875160e+05 );
    CHECK( params.b2   == -2.844167579e+04 );

    lines = {
        "    310.000    335.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 2.443570337D+02-1.399445548D+00 2.113509996D-03",
        " 0.000000000D+00 0.000000000D+00                -1.665032895D+05-1.059172219D+03" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  310.000 );
    CHECK( params.Tmax ==  335.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   ==  0.000000000e+00 );
    CHECK( params.a2   ==  0.000000000e+00 );
    CHECK( params.a3   ==  2.443570337e+02 );
    CHECK( params.a4   == -1.399445548e+00 );
    CHECK( params.a5   ==  2.113509996e-03 );
    CHECK( params.a6   ==  0.000000000e+00 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   == -1.665032895e+05 );
    CHECK( params.b2   == -1.059172219e+03 );

    lines = {
        "    335.000   2705.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-3.415474875D+05 0.000000000D+00 1.616932327D+01-1.517828471D-03 1.014852348D-06",
        " 0.000000000D+00 0.000000000D+00                -1.430478214D+05-8.374919535D+01" };

    params = parseNasaSpeciesThermoParams(lines);

    CHECK( params.Tmin ==  335.000 );
    CHECK( params.Tmax ==  2705.000 );
    CHECK( params.qN   ==  7 );
    CHECK( params.q1   == -2.0 );
    CHECK( params.q2   == -1.0 );
    CHECK( params.q3   ==  0.0 );
    CHECK( params.q4   ==  1.0 );
    CHECK( params.q5   ==  2.0 );
    CHECK( params.q6   ==  3.0 );
    CHECK( params.q7   ==  4.0 );
    CHECK( params.a1   == -3.415474875e+05 );
    CHECK( params.a2   ==  0.000000000e+00 );
    CHECK( params.a3   ==  1.616932327e+01 );
    CHECK( params.a4   == -1.517828471e-03 );
    CHECK( params.a5   ==  1.014852348e-06 );
    CHECK( params.a6   ==  0.000000000e+00 );
    CHECK( params.a7   ==  0.000000000e+00 );
    CHECK( params.b1   == -1.430478214e+05 );
    CHECK( params.b2   == -8.374919535e+01 );
}
