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

#include "ReactionEquation.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

ReactionEquation::ReactionEquation()
{}

ReactionEquation::ReactionEquation(Pairs<Species, double> const& species)
: m_species(species)
{}

ReactionEquation::ReactionEquation(const String& equation)
{
    // Split the reaction equation into two words: reactants and products
    auto two_words = split(equation, "=");

    // Assert the equation has a single equal sign `=`
    error(two_words.size() != 2,
        "Cannot parse the reaction equation `" +  equation + "`. ",
        "Expecting an equation with a single equal sign `=` separating "
        "reactants from products.");

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
        m_species.emplace_back(species, -number); // negative sign for reactants
    }

    // Iterave over all strings representing pair number and species name in the products
    for(auto word : products)
    {
        if(word == "+") continue;
        auto pair = split(word, "*");
        auto number = pair.size() == 2 ? tofloat(pair[0]) : 1.0;
        auto species = pair.size() == 2 ? pair[1] : pair[0];
        m_species.emplace_back(species, number); // positive sign for products
    }
}

ReactionEquation::ReactionEquation(const char* equation)
: ReactionEquation(String(equation))
{}

auto ReactionEquation::empty() const -> bool
{
    return m_species.empty();
}

auto ReactionEquation::size() const -> Index
{
    return m_species.size();
}

auto ReactionEquation::species() const -> Vec<Species>
{
    return vectorize(m_species, RKT_LAMBDA(x, x.first));
}

auto ReactionEquation::coefficients() const -> Vec<double>
{
    return vectorize(m_species, RKT_LAMBDA(x, x.second));
}

auto ReactionEquation::coefficient(const String& name) const -> double
{
    const auto idx = indexfn(m_species, RKT_LAMBDA(x, x.first.name() == name));
    return idx < size() ? m_species[idx].second : 0.0;
}

ReactionEquation::operator String() const
{
    auto coeffstr = [](double x) { return (x == 1.0) ? "" : std::to_string(x) + "*"; };
    String str;
    for(auto [x, y] : m_species) if(y < 0.0) str += coeffstr(-y) + x.name() + " + ";
    str = str.substr(0, str.size() - 3); // remove the leading string " + "
    str += "= ";
    for(auto [x, y] : m_species) if(y > 0.0) str += coeffstr(y) + x.name() + " ";
    str = str.substr(0, str.size() - 3); // remove the leading string " + "
    return str;
}

auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&
{
    out << String(equation);
    return out;
}

auto operator<(const ReactionEquation& lhs, const ReactionEquation& rhs) -> bool
{
    return lhs.species().front() < rhs.species().front();
}

auto operator==(const ReactionEquation& lhs, const ReactionEquation& rhs) -> bool
{
    // Check they have the same species
    if( !identical(lhs.species(), rhs.species()) )
        return false;

    // Check they have the same stoichiometric coefficients
    const auto eps = std::numeric_limits<double>::epsilon();
    for(auto s : lhs.species())
        if( std::abs(lhs.coefficient(s.name()) - rhs.coefficient(s.name())) > eps )
            return false;

    return true;
}

} // namespace Reaktoro
