// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "ReactionEquation.hpp"

// C++ includes
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
// Calcite + H+ = Ca+2 + HCO3-
// Dolomite + 2*H+ = Ca+2 + Mg+2 + 2*HCO3-
ReactionEquation::ReactionEquation()
{}

ReactionEquation::ReactionEquation(std::string equation)
: equation_str(equation)
{
    // Split the reaction equation into two words: reactants and products
    auto two_words = split(equation_str, "=");

    // Assert the equation has a single equal sign `=`
    Assert(two_words.size() == 2,
        "Cannot parse the reaction equation `" +  equation + "`.",
        "Expecting an equation with a single equal sign `=` separating "
        "reactants from products");

    // The reactants and products as string
    const auto& reactants_str = two_words[0];
    const auto& products_str = two_words[1];

    // Split the string representing the reactants and products at each `+` sign
    auto reactants = split(reactants_str, " ");
    auto products = split(products_str, " ");

    // Iterave over all strings representing pair number and species name in the reactants
    for(auto word : reactants)
    {
        if(word == "+") continue;
        auto pair = split(word, "*");
        auto number = pair.size() == 2 ? tofloat(pair[0]) : 1.0;
        auto species = pair.size() == 2 ? pair[1] : pair[0];
        equation_map.emplace(species, -number); // negative sign for reactants
    }

    // Iterave over all strings representing pair number and species name in the products
    for(auto word : products)
    {
        if(word == "+") continue;
        auto pair = split(word, "*");
        auto number = pair.size() == 2 ? tofloat(pair[0]) : 1.0;
        auto species = pair.size() == 2 ? pair[1] : pair[0];
        equation_map.emplace(species, number); // positive sign for products
    }
}

ReactionEquation::ReactionEquation(const std::map<std::string, double>& equation)
: equation_map(equation)
{
    std::stringstream reactants;
    std::stringstream products;
    for(auto pair : equation)
    {
        auto species = pair.first;
        auto stoichiometry = pair.second;
        if(stoichiometry == -1)
            reactants << species << " + ";
        else if(stoichiometry == +1)
            products << species << " + ";
        else if(stoichiometry < 0)
            reactants << -stoichiometry << "*" << species << " + ";
        else
            products << stoichiometry << "*" << species << " + ";
    }

    std::string reactants_str = reactants.str();
    std::string products_str = products.str();
    reactants_str.erase(reactants_str.end() - 3); // remove the leading " + "
    products_str.erase(products_str.end() - 3); // remove the leading " + "

    equation_str = reactants_str + " = " + products_str;
}

auto ReactionEquation::numSpecies() const -> unsigned
{
    return equation_map.size();
}

auto ReactionEquation::stoichiometry(std::string species) const -> double
{
    auto iter = equation_map.find(species);
    return iter != equation_map.end() ? iter->second : 0.0;
}

auto ReactionEquation::equation() const -> const std::map<std::string, double>&
{
    return equation_map;
}

ReactionEquation::operator std::string() const
{
    return equation_str;
}

auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&
{
    out << std::string(equation);
    return out;
}

} // namespace Reaktoro
