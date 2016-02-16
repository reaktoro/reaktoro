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
#include <set>
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

// Interpreter includes
#include "Keywords.hpp"
#include "Operators.hpp"
#include "Yaml.hpp"

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

/// Return a line with `Recipe` followed by `|` if needed.
auto fixRecipe(std::string line) -> std::string
{
    auto i = line.find("Recipe:");

    if(i == std::string::npos)
        return line; // there is no `Recipe:` in the line

    auto j = line.find('#');

    if(j < i)
        return line; // `Recipe:` appears after comment symbol #

    auto k = line.find_first_not_of(" ", i + 8);

    if(k < j) // `Recipe:` is followed by an inline list of compounds and their amounts
        return line;

    // Return modified line with `Recipe: |`
    return line.substr(0, i + 8) + " |" + line.substr(i + 8);
}

/// Return all entity names it can find in list of EntityValueUnits objects.
template<typename Triplet>
auto collectEntities(std::set<std::string>& list, const std::vector<Triplet>& triplets) -> void
{
    for(const Triplet& t : triplets)
        list.insert(t.entity);
}

/// Return all titrant names it can find in a list of EquilibriumConstraintNode objects.
template<typename Constraint>
auto collectTitrants(std::set<std::string>& list, const std::vector<Constraint>& constraints) -> void
{
    for(const Constraint& c : constraints)
    {
        list.insert(c.entity);
        list.insert(c.titrant1);
        list.insert(c.titrant2);
    }
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

auto collectCompoundsInEquilibriumNode(std::set<std::string>& set, const Node& node) -> void
{
    kwd::EquilibriumProblem kwd; node >> kwd;
    collectEntities(set, kwd.recipe);
    collectTitrants(set, kwd.ph);
    collectTitrants(set, kwd.species_amounts);
    collectTitrants(set, kwd.species_activities);
    collectTitrants(set, kwd.species_fugacities);
    collectEntities(set, kwd.inert_species);
    set.erase("");
}

auto collectCompoundsInSpeciationNode(std::set<std::string>& set, const Node& node) -> void
{
    kwd::SpeciationProblem kwd; node >> kwd;
    collectEntities(set, kwd.concentrations);
    collectTitrants(set, kwd.ph);
    collectTitrants(set, kwd.species_activities);
    collectTitrants(set, kwd.species_fugacities);
    collectEntities(set, kwd.inert_species);
    set.erase("");
}

auto collectCompoundsInAllEquilibriumNodes(std::set<std::string>& set, const Node& root) -> void
{
    for(auto child : root)
    {
        std::string key = lowercase(keyword(child));
        if(key == "equilibrium" || key == "equilibriumproblem")
            collectCompoundsInEquilibriumNode(set, child);
    }
}

auto collectCompoundsInAllSpeciationNodes(std::set<std::string>& set, const Node& root) -> void
{
    for(auto child : root)
    {
        std::string key = lowercase(keyword(child));
        if(key == "speciation" || key == "speciationproblem")
            collectCompoundsInSpeciationNode(set, child);
    }
}

} // namespace Reaktoro
