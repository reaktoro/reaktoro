// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Utils.hpp"

// C++ includes
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace {

// Return a line with added `- ` right before the first non-space char in the line.
auto fixHyphen(std::string line) -> std::string
{
    // Get the index of first char that is not a space
    auto i = line.find_first_not_of(" ");

    // Check if the line is a comment line, an empty line, or a line with spaces only
    if(i == std::string::npos || line[i] == '#')
        return line;

    // Add `- ` right before the first non-space char in the line
    return line.substr(0, i) + "- " + line.substr(i);
}

/// Return a line with `Mixture` followed by `|` if needed.
auto fixMixture(std::string line) -> std::string
{
    auto i = line.find("Mixture:");

    if(i == std::string::npos)
        return line; // there is no `Mixture:` in the line

    auto j = line.find('#');

    if(j < i)
        return line; // `Mixture:` appears after comment symbol #

    auto k = line.find_first_not_of(" ", i + 8);

    if(k < j) // `Mixture:` is followed by an inline mixture recipe
        return line;

    // Return modified line with `Mixture: |`
    return line.substr(0, i + 8) + " |" + line.substr(i + 8);
}

} // namespace

auto preprocess(std::string script) -> std::string
{
    std::istringstream iss(script);
    return preprocess(iss);
}

auto preprocess(std::istream& stream) -> std::string
{
    std::stringstream ss;
    std::string line;
    while(std::getline(stream, line))
    {
        line = fixHyphen(line);
        ss << line << std::endl;
    }
    return ss.str();
}

} // namespace Reaktoro
