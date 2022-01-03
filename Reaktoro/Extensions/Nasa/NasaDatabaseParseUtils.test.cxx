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

// C++ includes
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Extensions/Nasa/NasaDatabaseParseUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing NasaDatabaseParseUtils module", "[NasaDatabaseParseUtils]")
{
    using namespace NasaUtils;

    //======================================================================
    // Testing method NasaUtils::convertFortranScientificNumberToDouble
    //======================================================================

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

    //======================================================================
    // Testing method NasaUtils::getStringBetweenColumns
    //======================================================================

    const auto numbers = "8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06";

    CHECK( getStringBetweenColumns(numbers,  1, 15) == "8.906039290D+04"  );
    CHECK( getStringBetweenColumns(numbers, 16, 31) == "-9.750803930D+02" );
    CHECK( getStringBetweenColumns(numbers, 33, 33) == "5"                );

    //======================================================================
    // Testing method NasaUtils::getIntegerBetweenColumns
    //======================================================================

    CHECK( getIntegerBetweenColumns(numbers, 3, 4) == 90 );
    CHECK( getIntegerBetweenColumns(numbers, 1, 1) == 8  );

    //======================================================================
    // Testing method NasaUtils::getDoubleBetweenColumns
    //======================================================================

    CHECK( getDoubleBetweenColumns(numbers,  1, 15) == Approx(8.906039290e+04)  );
    CHECK( getDoubleBetweenColumns(numbers, 16, 31) == Approx(-9.750803930e+02) );
    CHECK( getDoubleBetweenColumns(numbers, 33, 33) == Approx(5.0)              );

    //======================================================================
    // Testing method NasaUtils::isCommentLine
    //======================================================================

    CHECK( isCommentLine("! Some text") );
    CHECK( isCommentLine("# Some text") );
    CHECK( isCommentLine("    ! Some text not trimmed ") );
    CHECK( isCommentLine("    # Some text not trimmed ") );

    //======================================================================
    // Testing method NasaUtils::parseFormula
    //======================================================================

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

    //======================================================================
    // Testing method NasaUtils::parseNasaSpeciesThermoParams
    //======================================================================

    Strings lines;
    NasaSpeciesThermoParams params;

    const Strings lines00 = {
        "    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10764.801",
        " 8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06",
        " 7.219318660D-10 7.329846740D-14                -3.425302600D+04-9.387539141D+00" };

    params = parseNasaSpeciesThermoParams(lines00);

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

    const Strings lines01 = {
        "    287.700    288.5007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128",
        "-1.494237648D+11 0.000000000D+00 8.652205900D+06-3.521604880D+04 3.968501633D+01",
        " 0.000000000D+00 0.000000000D+00                -1.866195945D+09-4.140143280D+07" };

    params = parseNasaSpeciesThermoParams(lines01);

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

    const Strings lines02 = {
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10772.157",
        "-3.834452340D+04-9.563261140D+02 8.200962780D+00-2.062130802D-04 9.872288990D-09",
        " 8.158366760D-12-7.527519660D-16                -3.423564020D+04-2.224772278D+01" };

    params = parseNasaSpeciesThermoParams(lines02);

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

    const Strings lines03 = {
        "   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428",
        " 2.291781634D+09-1.608862960D+06 4.312466360D+02-5.396508990D-02 3.531856210D-06",
        "-1.164403850D-10 1.527134223D-15                 1.258651434D+07-3.692101610D+03" };

    params = parseNasaSpeciesThermoParams(lines03);

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

    const Strings lines04 = {
        "    298.150    700.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 5.676714650D+06-6.561919460D+04 2.971622148D+02-6.132148300D-01 6.373067250D-04",
        "-2.569754132D-07 0.000000000D+00                 3.113488256D+05-1.716694497D+03" };

    params = parseNasaSpeciesThermoParams(lines04);

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

    const Strings lines05 = {
        "    700.000   1500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 8.246680410D+05-1.965549089D+03 5.659967790D+00 1.366241998D-02-1.241835551D-05",
        " 6.313550760D-09-1.192405400D-12                -3.594772610D+03-3.071970867D+01" };

    params = parseNasaSpeciesThermoParams(lines05);

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

    const Strings lines06 = {
        "   1500.000   2500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 4.325258320D+09-1.355254539D+07 1.757503878D+04-1.205560629D+01 4.623646210D-03",
        "-9.389890050D-07 7.891119490D-11                 8.506076540D+07-1.227417790D+05" };

    params = parseNasaSpeciesThermoParams(lines06);

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

    const Strings lines07 = {
        "    200.000    306.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-2.112979833D+06 0.000000000D+00 2.207431906D+02-1.164973854D+00 1.853971249D-03",
        " 0.000000000D+00 0.000000000D+00                -1.746831159D+05-9.949025210D+02" };

    params = parseNasaSpeciesThermoParams(lines07);

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

    const Strings lines08 = {
        "    306.000    310.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 6.705915562D+03-4.303760534D+01 6.919229155D-02",
        " 0.000000000D+00 0.000000000D+00                -8.349875160D+05-2.844167579D+04" };

    params = parseNasaSpeciesThermoParams(lines08);

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

    const Strings lines09 = {
        "    310.000    335.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        " 0.000000000D+00 0.000000000D+00 2.443570337D+02-1.399445548D+00 2.113509996D-03",
        " 0.000000000D+00 0.000000000D+00                -1.665032895D+05-1.059172219D+03" };

    params = parseNasaSpeciesThermoParams(lines09);

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

    const Strings lines10 = {
        "    335.000   2705.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        15300.000",
        "-3.415474875D+05 0.000000000D+00 1.616932327D+01-1.517828471D-03 1.014852348D-06",
        " 0.000000000D+00 0.000000000D+00                -1.430478214D+05-8.374919535D+01" };

    params = parseNasaSpeciesThermoParams(lines10);

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

    //======================================================================
    // Testing method NasaUtils::createNasaSpecies
    //======================================================================

    NasaSpecies species;

    const Strings lines11 = {
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
    };

    species = createNasaSpecies(lines11);

    CHECK( species.name == "MgO" );
    CHECK( species.comment == "Gurvich,1996a pt1 p398 pt2 p318." );
    CHECK( species.idcode == "tpis96" );
    CHECK( species.formula == Pairs<String, double>{{"MG", 1.0}, {"O", 1.0}} );
    CHECK( species.type == NasaSpeciesType::Gas );
    CHECK( species.molarmass == 40.3044000 );
    CHECK( species.dHf == 32261.307 );
    CHECK( species.dH0 == 8909.107 );
    CHECK( species.H0 == 0.0 );
    CHECK( species.Tmin == 200.0 );
    CHECK( species.Tmax == 20000.0 );
    CHECK( species.thermodata.size() == 3 );
    CHECK( species.thermodata[0] == parseNasaSpeciesThermoParams(StringsRange(lines11).segment(2, 3)) );
    CHECK( species.thermodata[1] == parseNasaSpeciesThermoParams(StringsRange(lines11).segment(5, 3)) );
    CHECK( species.thermodata[2] == parseNasaSpeciesThermoParams(StringsRange(lines11).segment(8, 3)) );

    const Strings lines12 = {
        "Mg(cr)            Hexagonal. Ref-Elm. Alcock,1993.",
        " 2 srd 93 MG  1.00    0.00    0.00    0.00    0.00 1   24.3050000          0.000",
        "    100.000    298.1507 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4979.161",
        "-5.412225134D+03 0.000000000D+00 1.458173723D+00 1.330204666D-02-4.098858502D-05",
        " 4.754339101D-08 0.000000000D+00                -7.759472010D+02-6.989702348D+00",
        "    298.150    923.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4979.161",
        "-2.860060304D+04 0.000000000D+00 3.398877384D+00-7.243962663D-04 1.405254188D-06",
        " 0.000000000D+00 0.000000000D+00                -1.089519906D+03-1.545973664D+01",
    };

    species = createNasaSpecies(lines12);

    CHECK( species.name == "Mg(cr)" );
    CHECK( species.comment == "Hexagonal. Ref-Elm. Alcock,1993." );
    CHECK( species.idcode == "srd 93" );
    CHECK( species.formula == Pairs<String, double>{{"MG", 1.0}} );
    CHECK( species.type == NasaSpeciesType::Condensed );
    CHECK( species.molarmass == 24.3050000 );
    CHECK( species.dHf == 0.000 );
    CHECK( species.dH0 == 4979.161 );
    CHECK( species.H0 == 0.0 );
    CHECK( species.Tmin == 100.0 );
    CHECK( species.Tmax == 923.0 );
    CHECK( species.thermodata.size() == 2 );
    CHECK( species.thermodata[0] == parseNasaSpeciesThermoParams(StringsRange(lines12).segment(2, 3)) );
    CHECK( species.thermodata[1] == parseNasaSpeciesThermoParams(StringsRange(lines12).segment(5, 3)) );

    const Strings lines13 = {
        "Jet-A(g)          McBride,1996. Faith,1971. Gracia-Salcedo,1988.          React.",
        " 2 g 8/01 C  12.00H  23.00    0.00    0.00    0.00 0  167.3110200    -249657.000",
        "    273.150   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        "-6.068695590D+05 8.328259590D+03-4.312321270D+01 2.572390455D-01-2.629316040D-04",
        " 1.644988940D-07-4.645335140D-11                -7.606962760D+04 2.794305937D+02",
        "   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0            0.000",
        " 1.858356102D+07-7.677219890D+04 1.419826133D+02-7.437524530D-03 5.856202550D-07",
        " 1.223955647D-11-3.149201922D-15                 4.221989520D+05-8.986061040D+02",
    };

    species = createNasaSpecies(lines13);

    CHECK( species.name == "Jet-A(g)" );
    CHECK( species.comment == "McBride,1996. Faith,1971. Gracia-Salcedo,1988.          React." );
    CHECK( species.idcode == "g 8/01" );
    CHECK( species.formula == Pairs<String, double>{{"C", 12.0}, {"H", 23.0}} );
    CHECK( species.type == NasaSpeciesType::Gas );
    CHECK( species.molarmass == 167.3110200 );
    CHECK( species.dHf == -249657.000 );
    CHECK( species.dH0 == 0.0 );
    CHECK( species.H0 == 0.0 );
    CHECK( species.Tmin == 273.150 );
    CHECK( species.Tmax == 6000.00 );
    CHECK( species.thermodata.size() == 2 );
    CHECK( species.thermodata[0] == parseNasaSpeciesThermoParams(StringsRange(lines13).segment(2, 3)) );
    CHECK( species.thermodata[1] == parseNasaSpeciesThermoParams(StringsRange(lines13).segment(5, 3)) );

    const Strings lines14 = {
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
    };

    species = createNasaSpecies(lines14);

    CHECK( species.name == "NaCN(II)" );
    CHECK( species.comment == "Lambda trans@288.5K. Chase,1998(3/66) pp631-3. Messer,1941." );
    CHECK( species.idcode == "g 8/01" );
    CHECK( species.formula == Pairs<String, double>{{"NA", 1.0}, {"C", 1.0}, {"N", 1.0}} );
    CHECK( species.type == NasaSpeciesType::Condensed );
    CHECK( species.molarmass == 49.0071700 );
    CHECK( species.dHf == -90709.000 );
    CHECK( species.dH0 == 19422.128 );
    CHECK( species.H0 == 0.0 );
    CHECK( species.Tmin == 197.700 );
    CHECK( species.Tmax == 288.500 );
    CHECK( species.thermodata.size() == 6 );
    CHECK( species.thermodata[0] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(2, 3)) );
    CHECK( species.thermodata[1] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(5, 3)) );
    CHECK( species.thermodata[2] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(8, 3)) );
    CHECK( species.thermodata[3] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(11, 3)) );
    CHECK( species.thermodata[4] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(14, 3)) );
    CHECK( species.thermodata[5] == parseNasaSpeciesThermoParams(StringsRange(lines14).segment(17, 3)) );

    const Strings lines15 = {
        "C2H2(L),acetyle   Acetylene. McBride,1996 pp84,92.",
        " 0 g 6/96 C   2.00H   2.00    0.00    0.00    0.00 1   26.0372800     207599.000",
        "    192.350      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
    };

    species = createNasaSpecies(lines15);

    CHECK( species.name == "C2H2(L),acetyle" );
    CHECK( species.comment == "Acetylene. McBride,1996 pp84,92." );
    CHECK( species.idcode == "g 6/96" );
    CHECK( species.formula == Pairs<String, double>{{"C", 2.0}, {"H", 2.0}} );
    CHECK( species.type == NasaSpeciesType::Condensed );
    CHECK( species.molarmass == 26.0372800 );
    CHECK( species.dHf == 0.0 );
    CHECK( species.dH0 == 0.0 );
    CHECK( species.H0 == 207599.000 );
    CHECK( species.Tmin == 192.350 );
    CHECK( species.Tmax == 192.350 );
    CHECK( species.thermodata.empty() );

    const Strings lines16 = {
        "N2O4(L)           Dinitrogen tetroxide. McBride,1996 pp85,93.",
        " 0 g 6/96 N   2.00O   4.00    0.00    0.00    0.00 1   92.0110000     -17549.000",
        "    298.150      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
    };

    species = createNasaSpecies(lines16);

    CHECK( species.name == "N2O4(L)" );
    CHECK( species.comment == "Dinitrogen tetroxide. McBride,1996 pp85,93." );
    CHECK( species.idcode == "g 6/96" );
    CHECK( species.formula == Pairs<String, double>{{"N", 2.0}, {"O", 4.0}} );
    CHECK( species.type == NasaSpeciesType::Condensed );
    CHECK( species.molarmass == 92.0110000 );
    CHECK( species.dHf == 0.0 );
    CHECK( species.dH0 == 0.0 );
    CHECK( species.H0 == -17549.000 );
    CHECK( species.Tmin == 298.150 );
    CHECK( species.Tmax == 298.150 );
    CHECK( species.thermodata.empty() );

    //======================================================================
    // Testing method NasaUtils::getNumberSpeciesBlocks
    //======================================================================
    Vec<NasaSpecies> speciesvec;

    const Strings lines17 = {
        "CLO3F(L)          Perchloryl Fluoride. McBride,1996 pp85,93.",
        " 0 g 6/96 CL  1.00O   3.00F   1.00    0.00    0.00 1  102.4496032     -47436.000",
        "    226.400      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "CL2(L)            Chlorine. McBride,1996 pp84,92.",
        " 0 g 6/96 CL  2.00    0.00    0.00    0.00    0.00 1   70.9060000     -22550.000",
        "    239.120      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "V2O5(L)           Liquid. Gurvich,1982 pt1 p68 pt2 p69.",
        " 1 tpis82 V   2.00O   5.00    0.00    0.00    0.00 2  181.8800000   -1551000.000",
        "    954.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        21210.000",
        " 0.000000000D+00 0.000000000D+00 2.285161723D+01 0.000000000D+00 0.000000000D+00",
        " 0.000000000D+00 0.000000000D+00                -1.880385445D+05-1.113109611D+02",
        "ZrC(cr)           Crystal. Chase,1998 pp658-60.",
        " 2 j12/64 ZR  1.00C   1.00    0.00    0.00    0.00 1  103.2347000    -196648.000",
        "    200.000    800.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        "-1.329596558D+05 1.766908163D+03-9.907099000D+00 5.427909170D-02-8.887722520D-05",
        " 7.332246190D-08-2.420411623D-11                -3.297124410D+04 5.279964170D+01",
        "    800.000   3805.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        " 8.297591180D+05-4.224196050D+03 1.308509211D+01-5.018095870D-03 2.155019459D-06",
        "-4.233046770D-10 3.266859200D-14                -7.112782350D+02-7.921945570D+01",
    };

    CHECK( getNumberSpeciesBlocks(lines14) == 1 );
    CHECK( getNumberSpeciesBlocks(lines15) == 1 );
    CHECK( getNumberSpeciesBlocks(lines16) == 1 );
    CHECK( getNumberSpeciesBlocks(lines17) == 4 );

    //======================================================================
    // Testing method NasaUtils::getNumberTextLinesForNextSpeciesBlock
    //======================================================================

    CHECK( getNumberTextLinesForNextSpeciesBlock(StringsRange(lines17).segment(0))  == 3 ); // from CLO3F(L) to CL2(L)
    CHECK( getNumberTextLinesForNextSpeciesBlock(StringsRange(lines17).segment(3))  == 3 ); // from CL2(L) to V2O5(L)
    CHECK( getNumberTextLinesForNextSpeciesBlock(StringsRange(lines17).segment(6))  == 5 ); // from V2O5(L) to ZrC(cr)
    CHECK( getNumberTextLinesForNextSpeciesBlock(StringsRange(lines17).segment(11)) == 8 ); // from ZrC(cr) to end

    //======================================================================
    // Testing method NasaUtils::createNasaSpeciesVector
    //======================================================================

    speciesvec = createNasaSpeciesVector(lines17);

    CHECK( speciesvec.size() == 4 );
    CHECK( speciesvec[0] == createNasaSpecies(StringsRange(lines17).segment(0, 3)) );
    CHECK( speciesvec[1] == createNasaSpecies(StringsRange(lines17).segment(3, 3)) );
    CHECK( speciesvec[2] == createNasaSpecies(StringsRange(lines17).segment(6, 5)) );
    CHECK( speciesvec[3] == createNasaSpecies(StringsRange(lines17).segment(11, 8)) );

    //======================================================================
    // Testing method NasaUtils::getTextLinesForProducts
    //======================================================================

    const Strings lines18 = {
        "!                                                                               ",
        "!                                                                               ",
        "thermo                                                                          ",
        "    200.00   1000.00   6000.00  20000.     9/09/04                              ",
        "CLO3F(L)          Perchloryl Fluoride. McBride,1996 pp85,93.                    ",
        " 0 g 6/96 CL  1.00O   3.00F   1.00    0.00    0.00 1  102.4496032     -47436.000",
        "    226.400      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "CL2(L)            Chlorine. McBride,1996 pp84,92.                               ",
        " 0 g 6/96 CL  2.00    0.00    0.00    0.00    0.00 1   70.9060000     -22550.000",
        "    239.120      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "END PRODUCTS                                                                    ",
        "V2O5(L)           Liquid. Gurvich,1982 pt1 p68 pt2 p69.                         ",
        " 1 tpis82 V   2.00O   5.00    0.00    0.00    0.00 2  181.8800000   -1551000.000",
        "    954.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        21210.000",
        " 0.000000000D+00 0.000000000D+00 2.285161723D+01 0.000000000D+00 0.000000000D+00",
        " 0.000000000D+00 0.000000000D+00                -1.880385445D+05-1.113109611D+02",
        "ZrC(cr)           Crystal. Chase,1998 pp658-60.                                 ",
        " 2 j12/64 ZR  1.00C   1.00    0.00    0.00    0.00 1  103.2347000    -196648.000",
        "    200.000    800.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        "-1.329596558D+05 1.766908163D+03-9.907099000D+00 5.427909170D-02-8.887722520D-05",
        " 7.332246190D-08-2.420411623D-11                -3.297124410D+04 5.279964170D+01",
        "    800.000   3805.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        " 8.297591180D+05-4.224196050D+03 1.308509211D+01-5.018095870D-03 2.155019459D-06",
        "-4.233046770D-10 3.266859200D-14                -7.112782350D+02-7.921945570D+01",
        "END REACTANTS                                                                   ",
    };

    CHECK( getTextLinesForProducts(lines18).size() == 6 );
    CHECK( getTextLinesForProducts(lines18)[0] == lines18[4] );
    CHECK( getTextLinesForProducts(lines18)[5] == lines18[9] );

    //======================================================================
    // Testing method NasaUtils::getTextLinesForReactants
    //======================================================================

    CHECK( getTextLinesForReactants(lines18).size() == 13 );
    CHECK( getTextLinesForReactants(lines18)[0]  == lines18[11] );
    CHECK( getTextLinesForReactants(lines18)[12] == lines18[23] );

    //======================================================================
    // Testing method NasaUtils::createTextLines
    //======================================================================

    const String text =
        "!                                                                               \n"
        "!                                                                               \n"
        "thermo                                                                          \n"
        "    200.00   1000.00   6000.00  20000.     9/09/04                              \n"
        "CLO3F(L)          Perchloryl Fluoride. McBride,1996 pp85,93.                    \n"
        " 0 g 6/96 CL  1.00O   3.00F   1.00    0.00    0.00 1  102.4496032     -47436.000\n"
        "    226.400      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000\n"
        "CL2(L)            Chlorine. McBride,1996 pp84,92.                               \n"
        "  # some out of place and nonsense comment preceded by spaces!                  \n"
        " 0 g 6/96 CL  2.00    0.00    0.00    0.00    0.00 1   70.9060000     -22550.000\n"
        "    239.120      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000\n"
        "END PRODUCTS                                                                    \n"
        "V2O5(L)           Liquid. Gurvich,1982 pt1 p68 pt2 p69.                         \n"
        " 1 tpis82 V   2.00O   5.00    0.00    0.00    0.00 2  181.8800000   -1551000.000\n"
        "    954.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        21210.000\n"
        " 0.000000000D+00 0.000000000D+00 2.285161723D+01 0.000000000D+00 0.000000000D+00\n"
        " 0.000000000D+00 0.000000000D+00                -1.880385445D+05-1.113109611D+02\n"
        "ZrC(cr)           Crystal. Chase,1998 pp658-60.                                 \n"
        " 2 j12/64 ZR  1.00C   1.00    0.00    0.00    0.00 1  103.2347000    -196648.000\n"
        "    200.000    800.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000\n"
        "-1.329596558D+05 1.766908163D+03-9.907099000D+00 5.427909170D-02-8.887722520D-05\n"
        "            ! Another silly comment and two empty lines below                   \n"
        "                                                                                \n"
        "\n"
        " 7.332246190D-08-2.420411623D-11                -3.297124410D+04 5.279964170D+01\n"
        "    800.000   3805.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000\n"
        " 8.297591180D+05-4.224196050D+03 1.308509211D+01-5.018095870D-03 2.155019459D-06\n"
        "-4.233046770D-10 3.266859200D-14                -7.112782350D+02-7.921945570D+01\n"
        "END REACTANTS                                                                   \n"
        "                                                                                \n";

    const Strings cleaned_and_right_trimmed_textlines = {
        "thermo",
        "    200.00   1000.00   6000.00  20000.     9/09/04",
        "CLO3F(L)          Perchloryl Fluoride. McBride,1996 pp85,93.",
        " 0 g 6/96 CL  1.00O   3.00F   1.00    0.00    0.00 1  102.4496032     -47436.000",
        "    226.400      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "CL2(L)            Chlorine. McBride,1996 pp84,92.",
        " 0 g 6/96 CL  2.00    0.00    0.00    0.00    0.00 1   70.9060000     -22550.000",
        "    239.120      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000",
        "END PRODUCTS",
        "V2O5(L)           Liquid. Gurvich,1982 pt1 p68 pt2 p69.",
        " 1 tpis82 V   2.00O   5.00    0.00    0.00    0.00 2  181.8800000   -1551000.000",
        "    954.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        21210.000",
        " 0.000000000D+00 0.000000000D+00 2.285161723D+01 0.000000000D+00 0.000000000D+00",
        " 0.000000000D+00 0.000000000D+00                -1.880385445D+05-1.113109611D+02",
        "ZrC(cr)           Crystal. Chase,1998 pp658-60.",
        " 2 j12/64 ZR  1.00C   1.00    0.00    0.00    0.00 1  103.2347000    -196648.000",
        "    200.000    800.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        "-1.329596558D+05 1.766908163D+03-9.907099000D+00 5.427909170D-02-8.887722520D-05",
        " 7.332246190D-08-2.420411623D-11                -3.297124410D+04 5.279964170D+01",
        "    800.000   3805.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         5862.000",
        " 8.297591180D+05-4.224196050D+03 1.308509211D+01-5.018095870D-03 2.155019459D-06",
        "-4.233046770D-10 3.266859200D-14                -7.112782350D+02-7.921945570D+01",
        "END REACTANTS",
    };

    std::stringstream ss;
    ss << text;

    CHECK( createTextLines(ss) == cleaned_and_right_trimmed_textlines);
}

