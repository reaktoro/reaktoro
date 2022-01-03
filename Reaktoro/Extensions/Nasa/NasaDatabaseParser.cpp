// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "NasaDatabaseParser.hpp"

// // C++ includes
// #include <iostream>

// // Reaktoro includes
// #include <Reaktoro/Common/Exception.hpp>
// #include <Reaktoro/Common/StringUtils.hpp>

// namespace Reaktoro {


// /// Convert a Fortran number in scientific format to `double`.
// auto convertFortranScientificNumberToDouble(String str) -> double
// {
//     const auto pos = str.find('D');
//     str.replace(pos, 1, 1, 'E');
//     return atof(str.c_str());
// }

// auto parseIntoLineStrings(std::istream& file) -> Strings
// {

// }

// /// Return true if a line in the Nasa database file is a comment line (it starts with ! or #).
// auto isCommentLine(const String& line) -> bool
// {
//     const auto firstchar = trimleft(line).front();
//     return firstchar == '!' || firstchar == '#';
// }

// /// Return the line position corresponding to a line that is not comment.
// auto skipCommentLines(const Strings& lines, Index linepos) -> Index
// {
//     const auto numlines = lines.size();
//     while(linepos < numlines && isCommentLine(lines[linepos]))
//         ++linepos;
//     return linepos;
// }

// auto createNextNasaSpecies(const Strings& lines, Index& linepos) -> NasaSpecies
// {
//     linepos = skipCommentLines(lines, linepos);

// }

// auto initNasaSpeciesBasicInfo(NasaSpecies& species, const Strings& lines, Index& linepos) -> void
// {
//     const auto record1 = lines[linepos];
//     species.name    = getStringBetweenColumns(record1,  1, 18);
//     species.comment = getStringBetweenColumns(record1, 19, 80);
//     linepos += 1;
//     const auto record2 = lines[linepos];
//     const auto numTintervals = getIntBetweenColumns(record2, 1, 2);

//     species.idcode = getStringBetweenColumns(record2, 4, 9);
// }

// /// Return the pairs of element symbols and coefficients in a Nasa-format formula string.
// auto parseFormula(const String& formula) -> Pairs<String, double>
// {
//     assert(formula.size() == 40);
//     Pairs<String, double> pairs;
//     auto offset = 0;
//     for(auto i = 0; i < 5; ++i)
//     {
//         const auto element = formula.substr(offset, 2);
//         const auto coeffstr = formula.substr(offset + 2, 6);
//         const auto coeff = tofloat(coeffstr);
//         if(coeff != 0.0)
//         {
//             pairs.push_back({ element, coeff });
//             error(!trim(element).empty(), "Cannot accept zero coefficient for element ", element, ".");
//         }
//     }
//     return pairs;
// }

// struct NasaDatabaseParser::Impl
// {

// };

// NasaDatabaseParser::NasaDatabaseParser()
// : pimpl(new Impl())
// {
// }

// NasaDatabaseParser::NasaDatabaseParser(const NasaDatabaseParser& other)
// : pimpl(new Impl(*other.pimpl))
// {

// }

// NasaDatabaseParser::~NasaDatabaseParser()
// {}

// auto NasaDatabaseParser::operator=(NasaDatabaseParser other) -> NasaDatabaseParser&
// {
//     pimpl = std::move(other.pimpl);
//     return *this;
// }


// } // namespace Reaktoro
