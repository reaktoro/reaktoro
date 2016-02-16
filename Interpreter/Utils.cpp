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

auto collectCompoundsInEquilibriumNode(std::set<std::string>& compounds, const Node& node) -> void
{
    kwd::EquilibriumProblem kwd; node >> kwd;
    collectEntities(compounds, kwd.recipe);
    collectTitrants(compounds, kwd.ph);
    collectTitrants(compounds, kwd.species_amounts);
    collectTitrants(compounds, kwd.species_activities);
    collectTitrants(compounds, kwd.species_fugacities);
    collectEntities(compounds, kwd.inert_species);
    compounds.erase("");
}

auto collectCompoundsInSpeciationNode(std::set<std::string>& compounds, const Node& node) -> void
{
    kwd::SpeciationProblem kwd; node >> kwd;
    collectEntities(compounds, kwd.concentrations);
    collectTitrants(compounds, kwd.ph);
    collectTitrants(compounds, kwd.species_activities);
    collectTitrants(compounds, kwd.species_fugacities);
    collectEntities(compounds, kwd.inert_species);
    compounds.erase("");
}

auto collectCompounds(const Node& root) -> std::vector<std::string>
{
    std::set<std::string> compounds;
    for(auto child : root)
    {
        std::string key = lowercase(keyword(child));
        if(key == "equilibrium" || key == "equilibriumproblem")
            collectCompoundsInEquilibriumNode(compounds, child);
        if(key == "speciation" || key == "speciationproblem")
            collectCompoundsInSpeciationNode(compounds, child);
    }
    return {compounds.begin(), compounds.end()};
}

auto collectElements(std::string compound, const Database& database, std::set<std::string>& elements) -> void
{
    if(database.containsAqueousSpecies(compound))
        for(auto e : database.aqueousSpecies(compound).elements()) // compound found as aqueous species
            elements.insert(e.first.name());
    else if(database.containsGaseousSpecies(compound))
        for(auto e : database.gaseousSpecies(compound).elements()) // compound found as gaseous species
            elements.insert(e.first.name());
    else if(database.containsMineralSpecies(compound))
        for(auto e : database.mineralSpecies(compound).elements()) // compound found as mineral species
            elements.insert(e.first.name());
    else for(auto e : Reaktoro::elements(compound)) // compound not found in the database, so it is a formula type
            elements.insert(e.first);
}

auto identifyElements(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::set<std::string> elements;
    for(auto compound : compounds)
        collectElements(compound, database, elements);
    return {elements.begin(), elements.end()};
}

auto filterGaseousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::set<std::string> gases;
    for(auto compound : compounds)
        if(database.containsGaseousSpecies(compound))
            gases.push_back(compound);
    return {gases.begin(), gases.end()};
}

auto filterMineralSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::set<std::string> minerals;
    for(auto compound : compounds)
        if(database.containsMineralSpecies(compound))
            minerals.push_back(compound);
    return {minerals.begin(), minerals.end()};
}

auto hasSpeciation(const Node& root) -> bool
{
    for(auto child : root)
    {
        std::string key = lowercase(keyword(child));
        if(key == "speciation" || key == "speciationproblem")
            return true;
    }
    return false;
}

auto speciate(const Node& root, const Database& database,
    std::vector<std::string>& aqueous_species,
    std::vector<std::string>& gaseous_species,
    std::vector<std::string>& mineral_species) -> void
{
    const auto has_speciation = hasSpeciation(root);
    const auto compounds = collectCompounds(root);
    const auto elements = identifyElements(compounds, database);

    aqueous_species = database.aqueousSpeciesWithElements(elements);
    gaseous_species = filterGaseousSpecies(compounds, database);
    mineral_species = filterMineralSpecies(compounds, database);

    if(hasSpeciation(root))
        mineral_species = database.mineralSpeciesWithElements(elements);
}

} // namespace Reaktoro
