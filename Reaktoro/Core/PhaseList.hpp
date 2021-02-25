// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

// Forward declaration of PhaseListBase
template<typename Data>
class PhaseListBase;

/// The specialized container to deal with a collection of Phase objects.
using PhaseList = PhaseListBase<Vec<Phase>>;

/// The specialized container to deal with a const reference view of a collection of Phase objects.
using PhaseListConstRef = PhaseListBase<const Vec<Phase>&>;

/// A type used as a collection of phases.
template<typename Data>
class PhaseListBase
{
public:
    /// Construct a default PhaseListBase object.
    PhaseListBase()
    {}

    /// Construct an PhaseListBase object with given phase.
    PhaseListBase(const Vec<Phase>& phase)
    : m_phases(phase)
    {}

    /// Construct an PhaseListBase object with given another one.
    template<typename OtherData>
    PhaseListBase(const PhaseListBase<OtherData>& other)
    : m_phases(other.m_phases)
    {}

    /// Append a new phase to the list of phase.
    auto append(const Phase& phase)
    {
        m_phases.push_back(phase);
    }

    /// Return the internal collection of Phase objects.
    auto data() const -> const Data&
    {
        return m_phases;
    }

    /// Return the number of phases in the collection.
    auto size() const -> Index
    {
        return m_phases.size();
    }

    /// Return the Phase object with given index.
    auto operator[](Index i) const -> const Phase&
    {
        return m_phases[i];
    }

    /// Return the index of the phase with given unique name or the number of phases if not found.
    auto find(const String& name) const -> Index
    {
        return findWithName(name);
    }

    /// Return the index of the phase with given unique name or the number of phases if not found.
    auto findWithName(const String& name) const -> Index
    {
        return indexfn(m_phases, RKT_LAMBDA(p, p.name() == name));
    }

    /// Return the index of the phase containing the species with given index or number of phases if not found.
    auto findWithSpecies(Index index) const -> Index
    {
        auto counter = 0;
        for(auto i = 0; i < size(); ++i) {
            counter += m_phases[i].species().size();
            if(counter > index)
                return i;
        }
        return size();
    }

    /// Return the index of the phase with given unique species name or the number of phases if not found.
    auto findWithSpecies(const String& name) const -> Index
    {
        return indexfn(m_phases, RKT_LAMBDA(p, containsfn(p.species(), RKT_LAMBDA(s, s.name() == name))));
    }

    /// Return the index of the first phase with given aggregate state of species or the number of phases if not found.
    auto findWithAggregateState(AggregateState option) const -> Index
    {
        return indexfn(m_phases, RKT_LAMBDA(p, p.aggregateState() == option));
    }

    /// Return the index of the first phase with given state of matter or the number of phases if not found.
    auto findWithStateOfMatter(StateOfMatter option) const -> Index
    {
        return indexfn(m_phases, RKT_LAMBDA(p, p.stateOfMatter() == option));
    }

    /// Return the index of the phase with given unique name or throw a runtime error if not found.
    auto index(const String& name) const -> Index
    {
        return indexWithName(name);
    }

    /// Return the index of the phase with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index
    {
        const auto idx = findWithName(name);
        error(idx >= size(), "Could not find any Phase object with name ", name, ".");
        return idx;
    }

    /// Return the index of the phase containing the species with given index or throw a runtime error if not found.
    auto indexWithSpecies(Index index) const -> Index
    {
        const auto idx = findWithSpecies(index);
        error(idx >= size(), "Could not find any Phase object containing a Species object with index ", index, ".");
        return idx;
    }

    /// Return the index of the phase with given unique species name or throw a runtime error if not found.
    auto indexWithSpecies(const String& name) const -> Index
    {
        const auto idx = findWithSpecies(name);
        error(idx >= size(), "Could not find any Phase object containing a Species object with name ", name, ".");
        return idx;
    }

    /// Return the index of the first phase with given aggregate state of species or throw a runtime error if not found.
    auto indexWithAggregateState(AggregateState option) const -> Index
    {
        const auto idx = findWithAggregateState(option);
        error(idx >= size(), "Could not find any Phase object whose species have aggregate state ", option, ".");
        return idx;
    }

    /// Return the index of the first phase with given state of matter or throw a runtime error if not found.
    auto indexWithStateOfMatter(StateOfMatter option) const -> Index
    {
        const auto idx = findWithStateOfMatter(option);
        error(idx >= size(), "Could not find any Phase object with state of matter ", option, ".");
        return idx;
    }

    /// Return all phases with given names.
    auto withNames(const StringList& names) const -> PhaseList
    {
        return vectorize(names, RKT_LAMBDA(name, m_phases[indexWithName(name)]));
    }

    /// Return all phases with given state of matter.
    auto withStateOfMatter(StateOfMatter option) const -> PhaseList
    {
        return filter(m_phases, RKT_LAMBDA(p, p.stateOfMatter() == option));
    }

    /// Return all phases whose species have the given aggregate state.
    auto withAggregateState(AggregateState option) const -> PhaseList
    {
        return filter(m_phases, RKT_LAMBDA(p, p.aggregateState() == option));
    }

    /// Return the number of species over all phases up to the one with given index.
    auto numSpeciesUntilPhase(Index iphase) const -> Index
    {
        auto sum = 0;
        for(auto i = 0; i < iphase; ++i)
            sum += m_phases[i].species().size();
        return sum;
    }

    /// Convert this PhaseListBase object into its Data.
    operator Data() { return m_phases; }

    /// Convert this PhaseListBase object into its Data.
    operator Data() const { return m_phases; }

private:
    /// The species stored in the list.
    Data m_phases;

public:
    /// Construct an PhaseList object with given begin and end iterators.
    template<typename InputIterator>
    PhaseListBase(InputIterator begin, InputIterator end) : m_phases(begin, end) {}

    /// Return begin const iterator of this PhaseList instance (for STL compatibility reasons).
    auto begin() const { return data().begin(); }

    /// Return begin iterator of this PhaseList instance (for STL compatibility reasons).
    auto begin() { return data().begin(); }

    /// Return end const iterator of this PhaseList instance (for STL compatibility reasons).
    auto end() const { return data().end(); }

    /// Return end iterator of this PhaseList instance (for STL compatibility reasons).
    auto end() { return data().end(); }

    /// Append a new Phase at the back of the container (for STL compatibility reasons).
    auto push_back(const Phase& species) -> void { append(species); }

    /// Insert a container of Phase objects into this PhaseList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_phases.insert(pos, begin, end); }

    /// The type of the value stored in a PhaseList (for STL compatibility reasons).
    using value_type = Phase;
};

} // namespace Reaktoro
