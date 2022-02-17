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

#include "PhaseList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>

namespace Reaktoro {

PhaseList::PhaseList()
{}

PhaseList::PhaseList(const Vec<Phase>& phase)
: m_phases(phase)
{}

auto PhaseList::append(const Phase& phase) -> void
{
    m_phases.push_back(phase);
}

auto PhaseList::data() const -> const Vec<Phase>&
{
    return m_phases;
}

auto PhaseList::empty() const -> bool
{
    return m_phases.empty();
}

auto PhaseList::size() const -> Index
{
    return m_phases.size();
}

auto PhaseList::species() const -> Vec<Species>
{
    auto num_species = 0;
    for(const auto& phase : m_phases)
        num_species += phase.species().size();
    Vec<Species> species;
    species.reserve(num_species);
    for(const auto& phase : m_phases)
        species = concatenate(species, phase.species().data());
    return species;
}

auto PhaseList::operator[](Index i) const -> const Phase&
{
    return m_phases[i];
}

auto PhaseList::operator[](Index i) -> Phase&
{
    return m_phases[i];
}

auto PhaseList::find(const String& name) const -> Index
{
    return findWithName(name);
}

auto PhaseList::findWithName(const String& name) const -> Index
{
    return indexfn(m_phases, RKT_LAMBDA(p, p.name() == name));
}

auto PhaseList::findWithSpecies(Index index) const -> Index
{
    auto counter = 0;
    for(auto i = 0; i < size(); ++i) {
        counter += m_phases[i].species().size();
        if(counter > index)
            return i;
    }
    return size();
}

auto PhaseList::findWithSpecies(const String& name) const -> Index
{
    return indexfn(m_phases, RKT_LAMBDA(p, containsfn(p.species(), RKT_LAMBDA(s, s.name() == name))));
}

auto PhaseList::findWithAggregateState(AggregateState option) const -> Index
{
    return indexfn(m_phases, RKT_LAMBDA(p, p.aggregateState() == option));
}

auto PhaseList::findWithStateOfMatter(StateOfMatter option) const -> Index
{
    return indexfn(m_phases, RKT_LAMBDA(p, p.stateOfMatter() == option));
}

auto PhaseList::index(const String& name) const -> Index
{
    return indexWithName(name);
}

auto PhaseList::indexWithName(const String& name) const -> Index
{
    const auto idx = findWithName(name);
    error(idx >= size(), "Could not find any Phase object with name ", name, ".");
    return idx;
}

auto PhaseList::indexWithSpecies(Index index) const -> Index
{
    const auto idx = findWithSpecies(index);
    error(idx >= size(), "Could not find any Phase object containing a Species object with index ", index, ".");
    return idx;
}

auto PhaseList::indexWithSpecies(const String& name) const -> Index
{
    const auto idx = findWithSpecies(name);
    error(idx >= size(), "Could not find any Phase object containing a Species object with name ", name, ".");
    return idx;
}

auto PhaseList::indexWithAggregateState(AggregateState option) const -> Index
{
    const auto idx = findWithAggregateState(option);
    error(idx >= size(), "Could not find any Phase object whose species have aggregate state ", option, ".");
    return idx;
}

auto PhaseList::indexWithStateOfMatter(StateOfMatter option) const -> Index
{
    const auto idx = findWithStateOfMatter(option);
    error(idx >= size(), "Could not find any Phase object with state of matter ", option, ".");
    return idx;
}

auto PhaseList::get(const String& name) const -> const Phase&
{
    return getWithName(name);
}

auto PhaseList::getWithName(const String& name) const -> const Phase&
{
    return m_phases[indexWithName(name)];
}

auto PhaseList::withNames(const StringList& names) const -> PhaseList
{
    return vectorize(names, RKT_LAMBDA(name, m_phases[indexWithName(name)]));
}

auto PhaseList::withStateOfMatter(StateOfMatter option) const -> PhaseList
{
    return filter(m_phases, RKT_LAMBDA(p, p.stateOfMatter() == option));
}

auto PhaseList::withAggregateState(AggregateState option) const -> PhaseList
{
    return filter(m_phases, RKT_LAMBDA(p, p.aggregateState() == option));
}

auto PhaseList::numSpeciesUntilPhase(Index iphase) const -> Index
{
    auto sum = 0;
    for(auto i = 0; i < iphase; ++i)
        sum += m_phases[i].species().size();
    return sum;
}

auto PhaseList::indicesPhasesArePure() const -> Indices
{
    Indices iphases;
    for(auto&& [i, phase] : enumerate(m_phases))
        if(phase.species().size() == 1)
            iphases.push_back(i);
    return iphases;
}

auto PhaseList::indicesPhasesAreSolution() const -> Indices
{
    Indices iphases;
    for(auto&& [i, phase] : enumerate(m_phases))
        if(phase.species().size() > 1)
            iphases.push_back(i);
    return iphases;
}

auto PhaseList::indicesSpeciesInPhases(const Indices& iphases) const -> Indices
{
    Indices res;
    for(auto i : iphases)
    {
        const auto ifirst = numSpeciesUntilPhase(i);
        const auto nspecies = m_phases[i].species().size();
        for(auto i = 0; i < nspecies; ++i)
            res.push_back(ifirst + i);
    }
    return res;
}

PhaseList::operator Vec<Phase>&()
{
    return m_phases;
}

PhaseList::operator Vec<Phase>const&() const
{
    return m_phases;
}

auto operator+(const PhaseList &a, const PhaseList &b) -> PhaseList
{
    return concatenate(a, b);
}

} // namespace Reaktoro
