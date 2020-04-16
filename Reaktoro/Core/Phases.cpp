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

#include "Phases.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

GenericPhase::GenericPhase()
{}

GenericPhase::GenericPhase(const StringList& species)
: names(species)
{}

GenericPhase::GenericPhase(const Speciate& elements)
: symbols(elements.symbols)
{}

GenericPhase::~GenericPhase()
{}

auto GenericPhase::setName(String name) -> GenericPhase&
{
    phasename = name;
    return *this;
}

auto GenericPhase::setStateOfMatter(StateOfMatter option) -> GenericPhase&
{
    stateofmatter = option;
    return *this;
}

auto GenericPhase::setAggregateState(AggregateState option) -> GenericPhase&
{
    aggregatestate = option;
    return *this;
}

auto GenericPhase::setActivityModel(const ActivityModel& model) -> GenericPhase&
{
    activitymodel = model;
    return *this;
}

auto GenericPhase::named(String name) -> GenericPhase&
{
    return setName(name);
}

auto GenericPhase::set(StateOfMatter option) -> GenericPhase&
{
    return setStateOfMatter(option);
}

auto GenericPhase::set(AggregateState option) -> GenericPhase&
{
    return setAggregateState(option);
}

auto GenericPhase::set(const ActivityModel& model) -> GenericPhase&
{
    return setActivityModel(model);
}

auto GenericPhase::name() const -> String
{
    return phasename;
}

auto GenericPhase::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GenericPhase::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GenericPhase::species() const -> const Strings&
{
    return names;
}

auto GenericPhase::elements() const -> const Strings&
{
    return symbols;
}

auto GenericPhase::activityModel() const -> const ActivityModel&
{
    return activitymodel;
}

auto GenericPhase::convert(const ThermoEngine& engine, const Strings& elements) const -> Phase
{
    error(aggregatestate == AggregateState::Undefined,
        "GenericPhase::convert requires an AggregateState value to be specified.\n"
        "Use method GenericPhase::setAggregateState to fix this.");

    auto species = engine.database().speciesWithAggregateState(aggregatestate);

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    Phase phase;
    phase = phase.withName(phasename);
    phase = phase.withSpecies(species);
    phase = phase.withStateOfMatter(stateofmatter);
    phase = phase.withStandardThermoPropsFn(engine.standardThermoPropsFn());
    phase = phase.withActivityPropsFn(activitymodel(species));

    return phase;
}


GenericPhases::GenericPhases()
{}

GenericPhases::GenericPhases(const StringList& species)
: names(species)
{}

GenericPhases::GenericPhases(const Speciate& elements)
: symbols(elements.symbols)
{}

GenericPhases::~GenericPhases()
{}

auto GenericPhases::setStateOfMatter(StateOfMatter option) -> GenericPhases&
{
    stateofmatter = option;
    return *this;
}

auto GenericPhases::setAggregateState(AggregateState option) -> GenericPhases&
{
    aggregatestate = option;
    return *this;
}

auto GenericPhases::setActivityModel(const ActivityModel& model) -> GenericPhases&
{
    activitymodel = model;
    return *this;
}

auto GenericPhases::set(StateOfMatter option) -> GenericPhases&
{
    return setStateOfMatter(option);
}

auto GenericPhases::set(AggregateState option) -> GenericPhases&
{
    return setAggregateState(option);
}

auto GenericPhases::set(const ActivityModel& model) -> GenericPhases&
{
    return setActivityModel(model);
}

auto GenericPhases::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GenericPhases::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GenericPhases::species() const -> const Strings&
{
    return names;
}

auto GenericPhases::elements() const -> const Strings&
{
    return symbols;
}

auto GenericPhases::activityModel() const -> const ActivityModel&
{
    return activitymodel;
}

auto GenericPhases::convert(const ThermoEngine& engine, const Strings& elements) const -> Vec<GenericPhase>
{
    error(aggregatestate != AggregateState::Undefined,
        "GenericPhases::convert requires an AggregateState value to be specified. "
        "Use method GenericPhases::set(AggregateState) to fix this.");

    auto species = engine.database().speciesWithAggregateState(aggregatestate);

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    Vec<GenericPhase> phases;
    phases.reserve(species.size());
    for(auto&& s : species)
    {
        GenericPhase phase(s.name());
        phase.setName(s.name());
        phase.setStateOfMatter(stateOfMatter());
        phase.setAggregateState(aggregateState());
        phase.setActivityModel(activityModel());

        phases.push_back( phase );
    }

    return phases;
}

} // namespace Reaktoro
