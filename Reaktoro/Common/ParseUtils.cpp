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

#include "ParseUtils.hpp"

// C++ includes
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

auto parseReaction(std::string reaction) -> std::map<std::string, double>
{
    std::map<std::string, double> equation;

    // Split the participating species in the reaction in words delimited by space
    std::vector<std::string> words = split(reaction, " ");

    // Create the pair entries
    for(const std::string& word : words)
    {
        std::vector<std::string> pair = split(word, ":");
        equation.emplace(pair[1], tofloat(pair[0]));
    }

    return equation;
}

} // namespace Reaktoro
