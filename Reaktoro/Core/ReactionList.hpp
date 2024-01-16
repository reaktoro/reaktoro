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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used as a collection of reactions.
class ReactionList
{
public:
    /// Construct a default ReactionList object.
    ReactionList();

    /// Construct an ReactionList object with given reactions.
    ReactionList(std::initializer_list<Reaction> reactions);

    /// Construct an ReactionList object with given reactions.
    ReactionList(const Vec<Reaction>& reactions);

    /// Append a new reaction to the list of reaction.
    auto append(const Reaction& reaction) -> void;

    /// Return the internal collection of Reaction objects.
    auto data() const -> const Vec<Reaction>&;

    /// Return true if there are no reactions in the collection.
    auto empty() const -> bool;

    /// Return the number of reactions in the collection.
    auto size() const -> Index;

    /// Return the Reaction object with given index.
    auto operator[](Index i) const -> const Reaction&;

    /// Return the Reaction object with given index.
    auto operator[](Index i) -> Reaction&;

    /// Return the index of the reaction with given unique name or the number of reactions if not found.
    auto find(const String& name) const -> Index;

    /// Return the index of the reaction with given unique name or the number of reactions if not found.
    auto findWithName(const String& name) const -> Index;

    /// Return the index of the reaction with given unique name or throw a runtime error if not found.
    auto index(const String& name) const -> Index;

    /// Return the index of the reaction with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index;

    /// Return the reaction with a given name.
    auto get(const String& name) const -> const Reaction&;

    /// Return the reaction with a given name.
    auto getWithName(const String& name) const -> const Reaction&;

    /// Return all reactions with given names.
    auto withNames(const StringList& names) const -> ReactionList;

    /// Convert this ReactionList object into its Vec<Reaction>.
    operator Vec<Reaction>&();

    /// Convert this ReactionList object into its Vec<Reaction>.
    operator Vec<Reaction>const&() const;

private:
    /// The reactions stored in the list.
    Vec<Reaction> m_reactions;

public:
    /// Construct an ReactionList object with given begin and end iterators.
    template<typename InputIterator>
    ReactionList(InputIterator begin, InputIterator end) : m_reactions(begin, end) {}

    /// Return begin const iterator of this ReactionList instance (for STL compatibility reasons).
    auto begin() const { return m_reactions.begin(); }

    /// Return begin iterator of this ReactionList instance (for STL compatibility reasons).
    auto begin() { return m_reactions.begin(); }

    /// Return end const iterator of this ReactionList instance (for STL compatibility reasons).
    auto end() const { return m_reactions.end(); }

    /// Return end iterator of this ReactionList instance (for STL compatibility reasons).
    auto end() { return m_reactions.end(); }

    /// Append a new Reaction at the back of the container (for STL compatibility reasons).
    auto push_back(const Reaction& species) -> void { append(species); }

    /// Insert a container of Reaction objects into this ReactionList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_reactions.insert(pos, begin, end); }

    /// The type of the value stored in a ReactionList (for STL compatibility reasons).
    using value_type = Reaction;
};

/// Return the concatenation of two ReactionList objects.
auto operator+(const ReactionList& a, const ReactionList& b) -> ReactionList;

} // namespace Reaktoro
