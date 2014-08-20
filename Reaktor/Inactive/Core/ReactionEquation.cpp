/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ReactionEquation.hpp"

// C++ includes
#include <sstream>

// Reaktor includes
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {

ReactionEquation::ReactionEquation()
{}

ReactionEquation::ReactionEquation(const std::string& equation)
{
    // Split the participating species in the reaction in words delimited by space
    std::vector<std::string> words = split(equation, " ");

    // Create the pair entries
    for(const std::string& word : words)
    {
        std::vector<std::string> pair = split(word, ":");
        push_back({pair[1], tofloat(pair[0])});
    }
}

ReactionEquation::ReactionEquation(const std::vector<std::string>& species, const std::vector<double>& stoichiometries)
{
    for(unsigned i = 0; i < species.size(); ++i)
        push_back({species[i], stoichiometries[i]});
}

ReactionEquation::operator std::string() const
{
    std::stringstream ss;

    for(const auto& pair : *this)
        ss << pair.second << ":" << pair.first << " ";

    return ss.str().substr(0, ss.str().size() - 1); // exclude the final characteres " "
}

auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&
{
    out << std::string(equation);

    return out;
}

} /* namespace Reaktor */
