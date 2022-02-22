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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used as a collection of phases.
class PhaseList
{
public:
    /// Construct a default PhaseList object.
    PhaseList();

    /// Construct an PhaseList object with given phase.
    PhaseList(const Vec<Phase>& phase);

    /// Append a new phase to the list of phase.
    auto append(const Phase& phase) -> void;

    /// Return the internal collection of Phase objects.
    auto data() const -> const Vec<Phase>&;

    /// Return true if there are no phases in the collection.
    auto empty() const -> bool;

    /// Return the number of phases in the collection.
    auto size() const -> Index;

    /// Return the species that compose the phases in the collection.
    auto species() const -> Vec<Species>;

    /// Return the Phase object with given index.
    auto operator[](Index i) const -> const Phase&;

    /// Return the Phase object with given index.
    auto operator[](Index i) -> Phase&;

    /// Return the index of the phase with given unique name or the number of phases if not found.
    auto find(const String& name) const -> Index;

    /// Return the index of the phase with given unique name or the number of phases if not found.
    auto findWithName(const String& name) const -> Index;

    /// Return the index of the phase containing the species with given index or number of phases if not found.
    auto findWithSpecies(Index index) const -> Index;

    /// Return the index of the phase with given unique species name or the number of phases if not found.
    auto findWithSpecies(const String& name) const -> Index;

    /// Return the index of the first phase with given aggregate state of species or the number of phases if not found.
    auto findWithAggregateState(AggregateState option) const -> Index;

    /// Return the index of the first phase with given state of matter or the number of phases if not found.
    auto findWithStateOfMatter(StateOfMatter option) const -> Index;

    /// Return the index of the phase with given unique name or throw a runtime error if not found.
    auto index(const String& name) const -> Index;

    /// Return the index of the phase with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index;

    /// Return the index of the phase containing the species with given index or throw a runtime error if not found.
    auto indexWithSpecies(Index index) const -> Index;

    /// Return the index of the phase with given unique species name or throw a runtime error if not found.
    auto indexWithSpecies(const String& name) const -> Index;

    /// Return the index of the first phase with given aggregate state of species or throw a runtime error if not found.
    auto indexWithAggregateState(AggregateState option) const -> Index;

    /// Return the index of the first phase with given state of matter or throw a runtime error if not found.
    auto indexWithStateOfMatter(StateOfMatter option) const -> Index;

    /// Return the phase with a given name.
    auto get(const String& name) const -> const Phase&;

    /// Return the phase with a given name.
    auto getWithName(const String& name) const -> const Phase&;

    /// Return all phases with given names.
    auto withNames(const StringList& names) const -> PhaseList;

    /// Return all phases with given state of matter.
    auto withStateOfMatter(StateOfMatter option) const -> PhaseList;

    /// Return all phases whose species have the given aggregate state.
    auto withAggregateState(AggregateState option) const -> PhaseList;

    /// Return the number of species over all phases up to the one with given index.
    auto numSpeciesUntilPhase(Index iphase) const -> Index;

    /// Return the indices of the phases with a single species.
    auto indicesPhasesArePure() const -> Indices;

    /// Return the indices of the phases with more than one species.
    auto indicesPhasesAreSolution() const -> Indices;

    /// Return the indices of the species in phases with a single species.
    auto indicesSpeciesInPurePhases() const -> Indices;

    /// Return the indices of the species in phases with more than one species.
    auto indicesSpeciesInSolutionPhases() const -> Indices;

    /// Return the indices of the species in the given phases.
    auto indicesSpeciesInPhases(const Indices& iphases) const -> Indices;

    /// Convert this PhaseList object into its Vec<Phase>.
    operator Vec<Phase>&();

    /// Convert this PhaseList object into its Vec<Phase>.
    operator Vec<Phase>const&() const;

private:
    /// The species stored in the list.
    Vec<Phase> m_phases;

public:
    /// Construct an PhaseList object with given begin and end iterators.
    template<typename InputIterator>
    PhaseList(InputIterator begin, InputIterator end) : m_phases(begin, end) {}

    /// Return begin const iterator of this PhaseList instance (for STL compatibility reasons).
    auto begin() const { return m_phases.begin(); }

    /// Return begin iterator of this PhaseList instance (for STL compatibility reasons).
    auto begin() { return m_phases.begin(); }

    /// Return end const iterator of this PhaseList instance (for STL compatibility reasons).
    auto end() const { return m_phases.end(); }

    /// Return end iterator of this PhaseList instance (for STL compatibility reasons).
    auto end() { return m_phases.end(); }

    /// Append a new Phase at the back of the container (for STL compatibility reasons).
    auto push_back(const Phase& species) -> void { append(species); }

    /// Insert a container of Phase objects into this PhaseList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_phases.insert(pos, begin, end); }

    /// The type of the value stored in a PhaseList (for STL compatibility reasons).
    using value_type = Phase;
};

/// Return the concatenation of two PhaseList objects.
auto operator+(const PhaseList &a, const PhaseList &b) -> PhaseList;

} // namespace Reaktoro
