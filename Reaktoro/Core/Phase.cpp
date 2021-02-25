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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace detail {

/// Return the molar masses of the species
auto molarMasses(const SpeciesList& species)
{
    ArrayXd molar_masses(species.size());
    transform(species, molar_masses, [](auto&& s) { return s.molarMass(); });
    return molar_masses;
}

/// Raise error if there is no common aggregate state for all species in the phase.
auto ensureCommonAggregateState(const SpeciesList& species)
{
    const auto aggregatestate = species[0].aggregateState();
    for(auto&& s : species)
        error(s.aggregateState() != aggregatestate,
            "The species in a phase need to have a common aggregate state.\n"
            "I got a list of species in which ", species[0].name(), " has\n"
            "aggregate state ", aggregatestate, " while ", s.name(), " has aggregate state ", s.aggregateState(), ".");
}

} // namespace detail

struct Phase::Impl
{
    /// The name of the phase
    String name;

    /// The state of matter of the phase.
    StateOfMatter state = StateOfMatter::Solid;

    /// The list of Species instances defining the phase
    SpeciesList species;

    /// The activity model function of the phase.
    ActivityPropsFn activity_props_fn;

    /// The ideal activity model function of the phase.
    ActivityPropsFn ideal_activity_props_fn;

    /// The molar masses of the species in the phase.
    ArrayXd species_molar_masses;
};

Phase::Phase()
: pimpl(new Impl())
{}

auto Phase::clone() const -> Phase
{
    Phase phase;
    *phase.pimpl = *pimpl;
    return phase;
}

auto Phase::withName(String name) -> Phase
{
    Phase copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Phase::withSpecies(const SpeciesList& species) -> Phase
{
    detail::ensureCommonAggregateState(species);
    Phase copy = clone();
    copy.pimpl->species = std::move(species);
    copy.pimpl->species_molar_masses = detail::molarMasses(copy.pimpl->species);
    return copy;
}

auto Phase::withStateOfMatter(StateOfMatter state) -> Phase
{
    Phase copy = clone();
    copy.pimpl->state = std::move(state);
    return copy;
}

auto Phase::withActivityPropsFn(ActivityPropsFn fn) -> Phase
{
    Phase copy = clone();
    copy.pimpl->activity_props_fn = std::move(fn);
    return copy;
}

auto Phase::withIdealActivityPropsFn(ActivityPropsFn fn) -> Phase
{
    Phase copy = clone();
    copy.pimpl->ideal_activity_props_fn = std::move(fn);
    return copy;
}

auto Phase::name() const -> String
{
    return pimpl->name;
}

auto Phase::stateOfMatter() const -> StateOfMatter
{
    return pimpl->state;
}

auto Phase::aggregateState() const -> AggregateState
{
    return species().size() ? species()[0].aggregateState() : AggregateState::Undefined;
}

auto Phase::species() const -> SpeciesListConstRef
{
    return pimpl->species;
}

auto Phase::species(Index idx) const -> const Species&
{
    return pimpl->species[idx];
}

auto Phase::activityPropsFn() const -> const ActivityPropsFn&
{
    return pimpl->activity_props_fn;
}

auto Phase::idealActivityPropsFn() const -> const ActivityPropsFn&
{
    return pimpl->ideal_activity_props_fn;
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
