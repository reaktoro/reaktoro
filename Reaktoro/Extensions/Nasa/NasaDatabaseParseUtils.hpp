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

/// The view to a sequence of strings in a vector of strings.
struct StringsRange
{
    using Iter = Strings::const_iterator;
    Iter b; ///< The begin iterator of the range.
    Iter e; ///< The end iterator of the range.
    auto operator[](Index i) const -> const String& { return *(b + i); }
    auto size() const { return e - b; }
    auto segment(Index pos, Index n) const { return StringsRange(b + pos, b + pos + n); }
    auto segment(Index pos) const { return segment(pos, size() - pos); }
    auto range(Index ibegin, Index iend) const { return StringsRange(b + ibegin, b + iend); }
    auto range(Index ibegin) const { return range(ibegin, size() - ibegin); }
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
/// This method will return a NasaThermoParams object initialized as follows:
///
/// ~~~c++
/// NasaThermoParams params;
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
auto parseNasaThermoParams(const StringsRange& lines) -> NasaThermoParams;

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
auto createNasaSpecies(const StringsRange& lines) -> NasaSpecies;

/// Create all chemical species with the given string lines comprising one or more species blocks.
/// This method assumes that `lines` only contain actual species data.
/// There are no comment lines (i.e., lines starting with `!` or `#`) and
/// lines starting with words `thermo`, `END PRODUCTS`, `END REACTANTS`.
auto createNasaSpeciesVector(const StringsRange& lines) -> Vec<NasaSpecies>;

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
