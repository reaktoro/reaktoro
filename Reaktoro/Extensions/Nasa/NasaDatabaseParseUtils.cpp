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

#include "NasaDatabaseParseUtils.hpp"

// C++ includes
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace NasaUtils {

auto convertFortranScientificNumberToDouble(const String& str) -> double
{
    const auto pos = str.find('D');
    if(pos == String::npos)
        return tofloat(str);
    auto copy = str;
    copy.replace(pos, 1, 1, 'e');
    return tofloat(copy);
}

auto getStringBetweenColumns(const String& str, Index begincol, Index endcol) -> String
{
    assert(begincol >= 1);
    assert(endcol >= begincol);
    return str.substr(begincol - 1, endcol - begincol + 1);
}

auto getIntegerBetweenColumns(const String& str, Index begincol, Index endcol) -> long
{
    const auto substr = getStringBetweenColumns(str, begincol, endcol);
    return std::stoi(substr);
}

auto getDoubleBetweenColumns(const String& str, Index begincol, Index endcol) -> double
{
    const auto substr = getStringBetweenColumns(str, begincol, endcol);
    return convertFortranScientificNumberToDouble(substr);
}

auto isCommentLine(const String& line) -> bool
{
    const auto pos = line.find_first_not_of(' ');
    const auto firstnonspacechar = line[pos];
    return firstnonspacechar == '!' || firstnonspacechar == '#';
}

auto parseFormula(const String& formula) -> Pairs<String, double>
{
    assert(formula.size() == 40);
    Pairs<String, double> pairs;
    auto offset = 0;
    for(auto i = 0; i < 5; ++i)
    {
        const auto element = trim(formula.substr(offset, 2));
        const auto coeffstr = formula.substr(offset + 2, 6);
        const auto coeff = tofloat(coeffstr);
        if(coeff != 0.0)
            pairs.push_back({ element, coeff });
        error(coeff == 0.0 && !element.empty(),
            "Cannot accept zero coefficient for element ", element, ".");
        offset += 8;
    }
    return pairs;
}

auto parseNasaThermoParams(const StringsRange& lines) -> NasaThermoParams
{
    assert(lines.size() == 3);

    const auto record3 = lines[0];
    const auto record4 = lines[1];
    const auto record5 = lines[2];

    const auto segment = split(getStringBetweenColumns(record3, 1, 22));
    error(segment.size() != 2, "Expecting only two double values, "
        "Tmin and Tmax, within the first 22 chars of line:\n", record3);

    const auto qvalues = split(getStringBetweenColumns(record3, 24, 63));
    error(qvalues.size() != 8, "Expecting 8 values for q exponents "
        "(last being zero) between columns 24 and 63 of line:\n", record3);

    NasaThermoParams params;
    params.Tmin = std::stod(segment[0]);
    params.Tmax = std::stod(segment[1]);
    params.qN = getIntegerBetweenColumns(record3, 23, 23);
    params.q1 = std::stod(qvalues[0]);
    params.q2 = std::stod(qvalues[1]);
    params.q3 = std::stod(qvalues[2]);
    params.q4 = std::stod(qvalues[3]);
    params.q5 = std::stod(qvalues[4]);
    params.q6 = std::stod(qvalues[5]);
    params.q7 = std::stod(qvalues[6]);
    params.a1 = getDoubleBetweenColumns(record4,  1, 16);
    params.a2 = getDoubleBetweenColumns(record4, 17, 32);
    params.a3 = getDoubleBetweenColumns(record4, 33, 48);
    params.a4 = getDoubleBetweenColumns(record4, 49, 64);
    params.a5 = getDoubleBetweenColumns(record4, 65, 80);
    params.a6 = getDoubleBetweenColumns(record5,  1, 16);
    params.a7 = getDoubleBetweenColumns(record5, 17, 32);
    params.b1 = getDoubleBetweenColumns(record5, 49, 64);
    params.b2 = getDoubleBetweenColumns(record5, 65, 80);

    error(params.qN != 7, "Cannot accept number of coefficients for "
        "Cp0 model different than 7. Got qN = ", params.qN, ".");

    return params;
}

auto convertIntegerToSpeciesType(Index value) -> NasaAggregateState
{
    return value == 0 ? NasaAggregateState::Gas : NasaAggregateState::Condensed;
}

auto getNumberTextLinesForNextSpeciesBlock(const StringsRange& lines) -> Index
{
    const auto record2 = lines[1];
    const auto numintervals = getIntegerBetweenColumns(record2, 1, 2);
    return numintervals == 0 ? 3 : 2 + 3 * numintervals;
}

auto getNumberSpeciesBlocks(const StringsRange& lines) -> Index
{
    auto counter = 0;
    for(const auto& line : lines)
        if(line.size() > 3 && line[0] == ' ' && std::isdigit(line[1]) && line[2] == ' ')
            ++counter;
    return counter;
}

auto createNasaSpecies(const StringsRange& lines) -> NasaSpecies
{
    assert(lines.size() >= 3);
    assert(lines.size() == getNumberTextLinesForNextSpeciesBlock(lines));

    NasaSpecies species;

    const auto record1 = lines[0];
    const auto record2 = lines[1];
    const auto record3 = lines[2];

    species.name      = trim(getStringBetweenColumns(record1,  1, 18));
    species.comment   = trim(getStringBetweenColumns(record1, 19, 80));
    species.idcode    = trim(getStringBetweenColumns(record2, 4, 9));
    species.formula   = parseFormula(getStringBetweenColumns(record2, 11, 50));
    species.aggregatestate      = convertIntegerToSpeciesType(getIntegerBetweenColumns(record2, 52, 52));
    species.molarmass = getDoubleBetweenColumns(record2, 53, 65);

    const auto numintervals = getIntegerBetweenColumns(record2, 1, 2);

    if(numintervals)
    {
        species.thermodata.resize(numintervals);
        for(auto i = 0; i < numintervals; ++i)
            species.thermodata[i] = parseNasaThermoParams(lines.segment(2 + 3*i, 3));

        species.dHf  = getDoubleBetweenColumns(record2, 66, 80);
        species.dH0  = getDoubleBetweenColumns(record3, 66, 80);
        species.H0   = 0.0;
        species.Tmin = species.thermodata.front().Tmin;
        species.Tmax = species.thermodata.back().Tmax;
    }
    else
    {
        species.dHf  = 0.0;
        species.dH0  = 0.0;
        species.H0   = getDoubleBetweenColumns(record2, 66, 80);
        species.Tmin = getDoubleBetweenColumns(record3, 1, 11);
        species.Tmax = species.Tmin;
    }

    return species;
}

auto createNasaSpeciesVector(const StringsRange& lines) -> Vec<NasaSpecies>
{
    Vec<NasaSpecies> vec;
    vec.reserve(getNumberSpeciesBlocks(lines));

    auto sublines = lines;
    while(sublines.size())
    {
        const auto numlines = getNumberTextLinesForNextSpeciesBlock(sublines);
        const auto nextlines = sublines.segment(0, numlines);
        vec.push_back(createNasaSpecies(nextlines));
        sublines = sublines.segment(numlines);
    }

    return vec;
}

auto getLineIndexBeginProducts(const StringsRange& lines) -> Index
{
    // This method goes from top to bottom and find the line containing the
    // word `thermo`. It then returns the index of this line incremented by 2,
    // since the first species block is after two lines. See below an example
    // in which species e- starts two lines after the `thermo` line:
    // ~~~
    // thermo
    //     200.00   1000.00   6000.00  20000.     9/09/04
    // e-                Ref-Species. Chase,1998 3/82.
    //  3 g12/98 E   1.00    0.00    0.00    0.00    0.00 0.000548579903          0.000
    //     298.150   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         6197.428
    //  0.000000000D+00 0.000000000D+00 2.500000000D+00 0.000000000D+00 0.000000000D+00
    //  0.000000000D+00 0.000000000D+00                -7.453750000D+02-1.172081224D+01
    // ~~~

    auto iline = 0;
    for(const auto& line : lines)
        if(line.substr(0, 6) == "thermo")
            break;
        else ++iline;
    iline += 2;
    return iline;
}

auto getLineIndexEndWord(const StringsRange& lines, const String& word) -> Index
{
    // This method goes from bottom to top (because there are fewer lines to
    // check for "END PRODUCTS", "END REACTANTS") to identify the line
    // containing the given word.

    const auto n = lines.size();
    auto i = 1;
    for(; i <= lines.size(); ++i)
        if(lines[n - i].substr(0, word.size()) == word)
            break;
    const auto iline = n - i;
    return iline;
}

auto getLineIndexEndProducts(const StringsRange& lines) -> Index
{
    return getLineIndexEndWord(lines, "END PRODUCTS");
}

auto getLineIndexBeginReactants(const StringsRange& lines) -> Index
{
    return getLineIndexEndProducts(lines) + 1; // next line after END PRODUCTS
}

auto getLineIndexEndReactants(const StringsRange& lines) -> Index
{
    return getLineIndexEndWord(lines, "END REACTANTS");
}

auto getTextLinesForProducts(const StringsRange& lines) -> StringsRange
{
    const auto ibegin = getLineIndexBeginProducts(lines);
    const auto iend = getLineIndexEndProducts(lines);

    error(ibegin >= lines.size(), "There is no `thermo` line in the database. "
        "Ensure a thermodynamic database file with NASA format has been provided.");

    error(iend >= lines.size(), "There is no `END PRODUCTS` line in the database. "
        "Ensure a thermodynamic database file with NASA format has been provided.");

    error(iend <= ibegin, "The `thermo` line is after the `END PRODUCTS` line. "
        "Ensure a thermodynamic database file with NASA format has been provided.");

    return lines.range(ibegin, iend);
}

auto getTextLinesForReactants(const StringsRange& lines) -> StringsRange
{
    const auto ibegin = getLineIndexBeginReactants(lines);
    const auto iend = getLineIndexEndReactants(lines);

    error(ibegin >= lines.size(), "There is no reactants in the database. "
        "Ensure a thermodynamic database file with NASA format has been provided.");

    error(iend >= lines.size(), "There is no `END REACTANTS` line in the database. "
        "Ensure a thermodynamic database file with NASA format has been provided.");

    return lines.range(ibegin, iend);
}

auto createTextLines(std::istream& file) -> Strings
{
    Strings lines;
    String line;
    while(std::getline(file, line))
    {
        line = trimright(line);
        if(line.empty() || isCommentLine(line))
            continue;
        lines.push_back(line);
    }
    return lines;
}

} // namespace NasaUtils
} // namespace Reaktoro
