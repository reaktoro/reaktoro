// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

ReactionEquation::ReactionEquation()
{}

ReactionEquation::ReactionEquation(Pairs<Species, double> const& species)
: m_species(species)
{}

ReactionEquation::ReactionEquation(const String& equation)
{
    for(const auto [name, coeff] : parseReactionEquation(equation))
        m_species.emplace_back(Species(name), coeff);
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
    auto coeffstr = [](double x) { return (x == 1.0) ? "" : str(x) + "*"; };
    String str;
    for(auto [x, y] : m_species) if(y < 0.0) str += coeffstr(-y) + x.name() + " + ";
    str = str.substr(0, str.size() - 3); // remove the leading string " + "
    str += " = ";
    for(auto [x, y] : m_species) if(y > 0.0) str += coeffstr(y) + x.name() + " + ";
    str = str.substr(0, str.size() - 3); // remove the leading string " + "
    if(str.substr(str.size() - 3, 3) == " = ") // check if leading string is " = ", which happens when the reaction has no explicit products (e.g., mineral-aqueous reaction when all aqueous species are in equilibrium)
        str = str.substr(0, str.size() - 3); // remove the leading string " = "
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
