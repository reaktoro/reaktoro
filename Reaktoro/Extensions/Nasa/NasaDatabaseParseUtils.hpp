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

#pragma once

// Reaktoro includes
#include <Reaktoro/Extensions/Nasa/NasaSpecies.hpp>

namespace Reaktoro {
namespace NasaUtils {

// Use this for testing.
// W(cr)             Crystal. Ref-Elm. Chase,1998 pp1925-8.
//  4 j 6/66 W   1.00    0.00    0.00    0.00    0.00 1  183.8400000          0.000
//     200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4973.000
// -6.824541400D+03-2.254249090D+02 4.976604610D+00-6.926436340D-03 1.202272986D-05
// -9.344133510D-09 2.818887123D-12                -3.510679270D+00-2.361334984D+01
//    1000.000   2600.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4973.000
//  5.530134840D+05-2.041485344D+03 5.870839470D+00-1.920714198D-03 1.067652983D-06
// -2.355109022D-10 2.160679310D-14                 1.163812518D+04-3.319171800D+01
//    2600.000   3200.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4973.000
//  2.474736879D+09 4.488921620D+06-1.235978300D+04 9.678565660D+00-3.556364610D-03
//  6.380420610D-07-4.521123450D-11                -2.029500909D+07 8.274369690D+04
//    3200.000   3680.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4973.000
// -1.755550399D+10 1.179059156D+07 1.177715365D+03-2.675166841D+00 7.252172480D-04
// -6.128007580D-08 0.000000000D+00                -9.702249190D+07-1.148926234D+03
// NaCN(II)          Lambda trans@288.5K. Chase,1998(3/66) pp631-3. Messer,1941.
//  6 g 8/01 NA  1.00C   1.00N   1.00    0.00    0.00 1   49.0071700     -90709.000
//     197.700    245.9007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
//  4.995073610D+08-1.007656780D+07 7.458386410D+04-2.134406974D+02-6.965856220D-02
//  1.540182072D-03-2.204709373D-06                 4.493712560D+07-3.984423000D+05
//     245.900    273.1007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
// -4.546907220D+08 4.339535660D+06 5.476070770D+02-1.324240852D+02 5.045446840D-01
// -5.791737950D-04 0.000000000D+00                -2.385188911D+07 3.106410058D+04
//     273.100    284.2007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
// -4.415096320D+08 0.000000000D+00 3.452631670D+04-1.661761169D+02 2.250608921D-01
//  0.000000000D+00 0.000000000D+00                -6.388306670D+06-1.596447272D+05
//     284.200    286.3007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
//  1.647379186D+09 0.000000000D+00-6.100370390D+04 1.429224260D+02 0.000000000D+00
//  0.000000000D+00 0.000000000D+00                 1.735057318D+07 3.142435378D+05
//     286.300    287.7007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
//  1.607424258D+11 0.000000000D+00-8.629992400D+06 3.290155880D+04-3.355904228D+01
//  0.000000000D+00 0.000000000D+00                 1.946284579D+09 4.175641230D+07
//     287.700    288.5007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        19422.128
// -1.494237648D+11 0.000000000D+00 8.652205900D+06-3.521604880D+04 3.968501633D+01
//  0.000000000D+00 0.000000000D+00                -1.866195945D+09-4.140143280D+07
// S(L)              Liquid. Ref-Elm. Gurvich,1989 pt1 p265 pt2 p160.
//  5 tpis89 S   1.00    0.00    0.00    0.00    0.00 3   32.0650000          0.000
//     388.360    428.1507 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4412.000
// -6.366550765D+07 0.000000000D+00 2.376860693D+03-7.888076026D+00 7.376076522D-03
//  0.000000000D+00 0.000000000D+00                -6.356594920D+05-1.186929589D+04
//     428.150    432.2507 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4412.000
//  0.000000000D+00 0.000000000D+00 6.928522306D+03-3.254655981D+01 3.824448176D-02
//  0.000000000D+00 0.000000000D+00                -9.832222680D+05-3.154806751D+04
//     432.250    453.1507 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4412.000
//  0.000000000D+00 0.000000000D+00 1.649945697D+02-6.843534977D-01 7.315907973D-04
//  0.000000000D+00 0.000000000D+00                -2.638846929D+04-7.681730097D+02
//     453.150    717.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4412.000
//  1.972984578D+06 0.000000000D+00-2.441009753D+01 6.090352889D-02-3.744069103D-05
//  0.000000000D+00 0.000000000D+00                 1.113013440D+04 1.363174183D+02
//     717.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         4412.000
//  0.000000000D+00 0.000000000D+00 3.848693429D+00 0.000000000D+00 0.000000000D+00
//  0.000000000D+00 0.000000000D+00                -8.284589830D+02-1.736128237D+01

/// The view to a sequence of strings in a vector of strings.
struct StringsRange
{
    using Iter = Strings::const_iterator;
    Iter b; ///< The begin iterator of the range.
    Iter e;   ///< The end iterator of the range.
    auto operator[](Index i) const -> const String& { return *(b + i); }
    auto size() const { return e - b; }
    auto range(Index pos, Index n) const { return StringsRange(b + pos, b + pos + n); }
    auto range(Index pos) const { return range(pos, size() - pos); }
    auto begin() const { return b; }
    auto end() const { return e; }
    StringsRange(const Iter& b, const Iter& e) : b(b), e(e) {}
    StringsRange(const Strings& vec) : b(vec.begin()), e(vec.end()) {}
};

/// Convert a Fortran number in scientific format to `double`.
auto convertFortranScientificNumberToDouble(const String& str) -> double;

/// Return the string between two given Fortran column numbers.
/// Example, `getStringBetweenColumns("AB CD EF", 4, 8)` returns `"CD EF"`.
auto getStringBetweenColumns(const String& str, Index begincol, Index endcol) -> String;

/// Return the integer value between two given Fortran column numbers.
/// Example, `getIntegerBetweenColumns("AB 972 EF", 4, 5)` returns `97`.
auto getIntegerBetweenColumns(const String& str, Index begincol, Index endcol) -> long;

/// Return the double value between two given Fortran column numbers.
/// Example, `getDoubleBetweenColumns("AB 9.72D+02 EF", 4, 11)` returns `9.72e+02`.
auto getDoubleBetweenColumns(const String& str, Index begincol, Index endcol) -> double;

/// Return true if a line in the Nasa database file is a comment line (it starts with ! or #).
auto isCommentLine(const String& line) -> bool;

/// Return the pairs of element symbols and coefficients in a Nasa-format formula string.
/// For example, species `ALH2F` has formula string:
///
/// `"AL  1.00H   2.00F   1.00    0.00    0.00"`.
///
/// For this formula string, the returned result of this method is:
///
/// `{{"AL", 1.0}, {"H", 2.0}, {"F", 1.0}}`.
auto parseFormula(const String& formula) -> Pairs<String, double>;

/// Return the Nasa thermodynamic paramerters in a block of lines.
/// Consider the following lines:
///
/// ~~~
///    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        10764.801
/// 8.906039290D+04-9.750803930D+02 5.870772320D+00 6.896663690D-03-4.045817280D-06
/// 7.219318660D-10 7.329846740D-14                -3.425302600D+04-9.387539141D+00
/// ~~~
///
/// This method will return a NasaSpeciesThermoParams object initialized as follows:
///
/// ~~~c++
/// NasaSpeciesThermoParams params;
/// params.Tmin =  200.0;
/// params.Tmax =  1000.0;
/// params.qN   =  7;
/// params.q1   = -2.0;
/// params.q2   = -1.0;
/// params.q3   =  0.0;
/// params.q4   =  1.0;
/// params.q5   =  2.0;
/// params.q6   =  3.0;
/// params.q7   =  4.0;
/// params.a1   =  8.906039290e+04;
/// params.a2   = -9.750803930e+02;
/// params.a3   =  5.870772320e+00;
/// params.a4   =  6.896663690e-03;
/// params.a5   = -4.045817280e-06;
/// params.a6   =  7.219318660e-10;
/// params.a7   =  7.329846740e-14;
/// params.b1   = -3.425302600e+04;
/// params.b2   = -9.387539141e+00;
/// ~~~
///
/// Note the values 10764.801 above was not used. This value is set
/// in @ref NasaSpecies::dH0 by @ref createNextSpecies when creating
/// a NasaSpecies object.
auto parseNasaSpeciesThermoParams(const StringsRange& lines) -> NasaSpeciesThermoParams;

/// Return the number of species blocks in a range of string lines.
/// This method assumes that the number of species is the number of lines
/// containing digits representing number of temperature invervals.
/// For example, for species Ag below, the line starting with " 3 "
/// represents a species block.
/// ----------------------------------------------------------------------------------
/// Ag                Hf:Cox,1989. Moore,1971. Gordon,1999.
///  3 g10/97 AG  1.00    0.00    0.00    0.00    0.00 0  107.8682000     284900.000
///     200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
auto getNumberSpeciesBlocks(const StringsRange& lines) -> Index;

/// Return the number of string lines in the next species block.
/// For example, if `lines` is:
/// ~~~
/// ALBr2             Gurvich,1996a pt1 p186 pt2 p149.
///  2 tpis96 AL  1.00BR  2.00    0.00    0.00    0.00 0  186.7895380    -140662.125
///     200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
///  3.199375870D+04-7.119178970D+02 9.478258110D+00-4.875531670D-03 5.516512990D-06
/// -3.340053040D-09 8.368476840D-13                -1.540591306D+04-1.742171366D+01
///    1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        13397.875
/// -3.523782900D+05 4.671544170D+02 7.111908190D+00-5.551709200D-04 3.166301130D-07
/// -5.521028330D-11 3.176725950D-15                -2.265004078D+04-2.695610360D+00
/// CH4(L)            Methane. McBride,1996 pp85,93.
///  0 g 6/96 C   1.00H   4.00    0.00    0.00    0.00 1   16.0424600     -89233.000
///     111.643      0.0000  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0            0.000
/// ~~~
/// then this method returns 8, i.e., the number of lines between `ALBr2` and `CH4(L)`.
auto getNumberTextLinesForNextSpeciesBlock(const StringsRange& lines) -> Index;

/// Create a chemical species with the given string lines.
auto createSpecies(const StringsRange& lines) -> NasaSpecies;

/// Create all chemical species with the given string lines comprising one or more species blocks.
/// This method assumes that `lines` only contain actual species data.
/// There are no comment lines (i.e., lines starting with `!` or `#`) and
/// lines starting with words `thermo`, `END PRODUCTS`, `END REACTANTS`.
auto createSpeciesVector(const StringsRange& lines) -> Vec<NasaSpecies>;

/// Return a range of text lines comprising the text blocks for product species.
/// This method returns all lines shown below:
/// ~~~
/// thermo
///     200.00   1000.00   6000.00  20000.     9/09/04
/// ...
/// ... the lines that go here are returned!
/// ...
/// END PRODUCTS
/// ~~~
auto getTextLinesForProducts(const StringsRange& lines) -> StringsRange;

/// Return a range of text lines comprising the text blocks for reactant species.
/// This method returns all lines shown below:
/// ~~~
/// END PRODUCTS
/// ...
/// ... the lines that go here are returned!
/// ...
/// END REACTANTS
/// ~~~
auto getTextLinesForReactants(const StringsRange& lines) -> StringsRange;

/// Return a vector with string lines comprising the contents of a database file.
/// In the process, all empty lines and comment lines (i.e., those starting
/// with `!` or `#`) are removed.
auto createTextLines(std::istream& file) -> Strings;

} // namemespace NasaUtils
} // namespace Reaktoro
