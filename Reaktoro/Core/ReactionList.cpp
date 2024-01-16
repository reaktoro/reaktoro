// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "ReactionList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>

namespace Reaktoro {

ReactionList::ReactionList()
{}

ReactionList::ReactionList(std::initializer_list<Reaction> reactions)
: m_reactions(std::move(reactions))
{}

ReactionList::ReactionList(const Vec<Reaction>& reactions)
: m_reactions(reactions)
{}

auto ReactionList::append(const Reaction& reaction) -> void
{
    m_reactions.push_back(reaction);
}

auto ReactionList::data() const -> const Vec<Reaction>&
{
    return m_reactions;
}

auto ReactionList::empty() const -> bool
{
    return m_reactions.empty();
}

auto ReactionList::size() const -> Index
{
    return m_reactions.size();
}

auto ReactionList::operator[](Index i) const -> const Reaction&
{
    return m_reactions[i];
}

auto ReactionList::operator[](Index i) -> Reaction&
{
    return m_reactions[i];
}

auto ReactionList::find(const String& name) const -> Index
{
    return findWithName(name);
}

auto ReactionList::findWithName(const String& name) const -> Index
{
    return indexfn(m_reactions, RKT_LAMBDA(p, p.name() == name));
}

auto ReactionList::index(const String& name) const -> Index
{
    return indexWithName(name);
}

auto ReactionList::indexWithName(const String& name) const -> Index
{
    const auto idx = findWithName(name);
    error(idx >= size(), "Could not find any Reaction object with name ", name, ".");
    return idx;
}

auto ReactionList::get(const String& name) const -> const Reaction&
{
    return getWithName(name);
}

auto ReactionList::getWithName(const String& name) const -> const Reaction&
{
    return m_reactions[indexWithName(name)];
}

auto ReactionList::withNames(const StringList& names) const -> ReactionList
{
    return vectorize(names, RKT_LAMBDA(name, m_reactions[indexWithName(name)]));
}

ReactionList::operator Vec<Reaction>&()
{
    return m_reactions;
}

ReactionList::operator Vec<Reaction>const&() const
{
    return m_reactions;
}

auto operator+(const ReactionList& a, const ReactionList& b) -> ReactionList
{
    return concatenate(a, b);
}

} // namespace Reaktoro
