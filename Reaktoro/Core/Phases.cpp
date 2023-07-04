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

#include "Phases.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>

namespace Reaktoro {

auto speciate(StringList const& substances) -> Speciate
{
    errorif(substances.empty(), "Expecting a non-empty list of substance formulas in method `speciate`.");
    Set<String> elements;
    for(auto const& substance : substances)
        for(auto const& [symbol, coeff] : parseChemicalFormula(substance))
            elements.insert(symbol);
    return Speciate{Strings(elements.begin(), elements.end())};
}

auto exclude(const StringList& tags) -> Exclude
{
    errorif(tags.empty(), "Expecting a non-empty list of tags in method `exclude`.");
    return Exclude{tags};
}

GeneralPhase::GeneralPhase()
{}

GeneralPhase::GeneralPhase(const StringList& species)
: names(species)
{}

GeneralPhase::GeneralPhase(const Speciate& elements)
: symbols(elements.symbols)
{}

GeneralPhase::GeneralPhase(const Speciate& elements, const Exclude& withtags)
: symbols(elements.symbols), excludetags(withtags.tags)
{}

GeneralPhase::GeneralPhase(const Exclude& withtags)
: excludetags(withtags.tags)
{}

GeneralPhase::~GeneralPhase()
{}

auto GeneralPhase::setName(String name) -> GeneralPhase&
{
    phasename = name;
    return *this;
}

auto GeneralPhase::setStateOfMatter(StateOfMatter option) -> GeneralPhase&
{
    stateofmatter = option;
    return *this;
}

auto GeneralPhase::setAggregateState(AggregateState option) -> GeneralPhase&
{
    aggregatestate = option;
    return *this;
}

auto GeneralPhase::setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GeneralPhase&
{
    other_aggregate_states = options;
    return *this;
}

auto GeneralPhase::setActivityModel(const ActivityModelGenerator& model) -> GeneralPhase&
{
    activity_model = model;
    return *this;
}

auto GeneralPhase::setIdealActivityModel(const ActivityModelGenerator& model) -> GeneralPhase&
{
    ideal_activity_model = model;
    return *this;
}

auto GeneralPhase::named(String name) -> GeneralPhase&
{
    return setName(name);
}

auto GeneralPhase::set(StateOfMatter option) -> GeneralPhase&
{
    return setStateOfMatter(option);
}

auto GeneralPhase::set(AggregateState option) -> GeneralPhase&
{
    return setAggregateState(option);
}

auto GeneralPhase::set(const ActivityModelGenerator& model) -> GeneralPhase&
{
    return setActivityModel(model);
}

auto GeneralPhase::name() const -> String
{
    return phasename;
}

auto GeneralPhase::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GeneralPhase::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GeneralPhase::additionalAggregateStates() const -> const Vec<AggregateState>&
{
    return other_aggregate_states;
}

auto GeneralPhase::species() const -> const Strings&
{
    return names;
}

auto GeneralPhase::elements() const -> const Strings&
{
    return symbols;
}

auto GeneralPhase::activityModel() const -> const ActivityModelGenerator&
{
    return activity_model;
}

auto GeneralPhase::idealActivityModel() const -> const ActivityModelGenerator&
{
    return ideal_activity_model;
}

auto GeneralPhase::convert(const Database& db, const Strings& elements) const -> Phase
{
    error(aggregatestate == AggregateState::Undefined,
        "GeneralPhase::convert requires an AggregateState value to be specified.\n"
        "Use method GeneralPhase::setAggregateState to fix this.");

    auto species = db.speciesWithAggregateState(aggregatestate);

    // If additional aggregate states provided, consider also other species in the database
    for(auto other_aggregate_state : other_aggregate_states)
    {
        auto other_species = db.speciesWithAggregateState(other_aggregate_state);
        if(other_species.size())
            species = concatenate(species, other_species);
    }

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    // Filter out species with provided tags in the exclude function
    if(excludetags.size())
        species = species.withoutTags(excludetags);

    errorif(species.empty(), "Expecting at least one species when defining a phase, but none was provided. Make sure you have listed the species names yourself or used the `speciate` method appropriately.")

    Phase phase;
    phase = phase.withName(phasename);
    phase = phase.withStateOfMatter(stateofmatter);
    phase = phase.withSpecies(species);
    phase = phase.withActivityModel(activity_model(species));
    phase = phase.withIdealActivityModel(ideal_activity_model(species));

    return phase;
}


GeneralPhasesGenerator::GeneralPhasesGenerator()
{}

GeneralPhasesGenerator::GeneralPhasesGenerator(const StringList& species)
: names(species)
{}

GeneralPhasesGenerator::GeneralPhasesGenerator(const Speciate& elements)
: symbols(elements.symbols)
{}

GeneralPhasesGenerator::GeneralPhasesGenerator(const Speciate& elements, const Exclude& withtags)
: symbols(elements.symbols), excludetags(withtags.tags)
{}

GeneralPhasesGenerator::GeneralPhasesGenerator(const Exclude& withtags)
: excludetags(withtags.tags)
{}

GeneralPhasesGenerator::~GeneralPhasesGenerator()
{}

auto GeneralPhasesGenerator::setStateOfMatter(StateOfMatter option) -> GeneralPhasesGenerator&
{
    stateofmatter = option;
    return *this;
}

auto GeneralPhasesGenerator::setAggregateState(AggregateState option) -> GeneralPhasesGenerator&
{
    aggregatestate = option;
    return *this;
}

auto GeneralPhasesGenerator::setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GeneralPhasesGenerator&
{
    other_aggregate_states = options;
    return *this;
}

auto GeneralPhasesGenerator::setActivityModel(const ActivityModelGenerator& model) -> GeneralPhasesGenerator&
{
    activity_model = model;
    return *this;
}

auto GeneralPhasesGenerator::setIdealActivityModel(const ActivityModelGenerator& model) -> GeneralPhasesGenerator&
{
    ideal_activity_model = model;
    return *this;
}

auto GeneralPhasesGenerator::set(StateOfMatter option) -> GeneralPhasesGenerator&
{
    return setStateOfMatter(option);
}

auto GeneralPhasesGenerator::set(AggregateState option) -> GeneralPhasesGenerator&
{
    return setAggregateState(option);
}

auto GeneralPhasesGenerator::set(const ActivityModelGenerator& model) -> GeneralPhasesGenerator&
{
    return setActivityModel(model);
}

auto GeneralPhasesGenerator::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GeneralPhasesGenerator::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GeneralPhasesGenerator::additionalAggregateStates() const -> const Vec<AggregateState>&
{
    return other_aggregate_states;
}

auto GeneralPhasesGenerator::species() const -> const Strings&
{
    return names;
}

auto GeneralPhasesGenerator::elements() const -> const Strings&
{
    return symbols;
}

auto GeneralPhasesGenerator::activityModel() const -> const ActivityModelGenerator&
{
    return activity_model;
}

auto GeneralPhasesGenerator::idealActivityModel() const -> const ActivityModelGenerator&
{
    return ideal_activity_model;
}

auto GeneralPhasesGenerator::convert(const Database& db, const Strings& elements) const -> Vec<GeneralPhase>
{
    error(aggregatestate == AggregateState::Undefined,
        "GeneralPhasesGenerator::convert requires an AggregateState value to be specified. "
        "Use method GeneralPhasesGenerator::set(AggregateState) to fix this.");

    auto species = db.speciesWithAggregateState(aggregatestate);

    // If additional aggregate states provided, consider also other species in the database
    for(auto other_aggregate_state : other_aggregate_states)
    {
        auto other_species = db.speciesWithAggregateState(other_aggregate_state);
        if(other_species.size())
            species = concatenate(species, other_species);
    }

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    // Filter out species with provided tags in the exclude function
    if(excludetags.size())
        species = species.withoutTags(excludetags);

    errorif(species.empty(), "Expecting at least one species when defining a list of single-species phases, but none was provided. Make sure you have listed the species names yourself or used the `speciate` method appropriately.")

    Vec<GeneralPhase> phases;
    phases.reserve(species.size());
    for(auto&& s : species)
    {
        GeneralPhase phase(s.name());
        phase.setName(s.name());
        phase.setStateOfMatter(stateOfMatter());
        phase.setAggregateState(aggregateState());
        phase.setAdditionalAggregateStates(additionalAggregateStates());
        phase.setActivityModel(activityModel());
        phase.setIdealActivityModel(idealActivityModel());

        phases.push_back( phase );
    }

    return phases;
}

Phases::Phases(const Database& db)
: db(db)
{}

auto Phases::add(const GeneralPhase& phase) -> void
{
    generalphases.push_back(phase);
}

auto Phases::add(const GeneralPhasesGenerator& generator) -> void
{
    generators.push_back(generator);
}

auto Phases::database() const -> const Database&
{
    return db;
}

auto Phases::generalPhases() const -> Vec<GeneralPhase> const&
{
    return generalphases;
}

auto Phases::generalPhasesGenerators() const -> Vec<GeneralPhasesGenerator> const&
{
    return generators;
}

auto Phases::convert() const -> Vec<Phase>
{
    // Return the element symbols in all stored GeneralPhase and GeneralPhasesGenerator objects.
    auto collect_all_element_symbols = [&]() -> Strings
    {
        Strings symbols;

        const auto collect_element_symbols_in_generalphase_or_generator = [&](const auto& phase)
        {
            Strings result;
            if(phase.elements().size())
                result = merge(result, phase.elements());
            if(phase.species().size())
                for(auto&& s : db.species().withNames(phase.species()))
                    result = merge(result, s.elements().symbols());
            if(phase.aggregateState() == AggregateState::Aqueous)
                result = merge(result, Strings{"H", "O"}); // ensure both H and O are considered in case there is aqueous phases
            return result;
        };

        for(const auto& phase : generalphases)
            symbols = merge(symbols, collect_element_symbols_in_generalphase_or_generator(phase));

        for(const auto& generator : generators)
            symbols = merge(symbols, collect_element_symbols_in_generalphase_or_generator(generator));

        return symbols;
    };

    // Return all given GeneralPhase objects together with the generated ones.
    auto collect_all_general_phases = [&](const Strings& symbols) -> Vec<GeneralPhase>
    {
        Vec<GeneralPhase> collected(generalphases);

        for(auto generator : generators)
        {
            const Vec<GeneralPhase> generated = generator.convert(db, symbols);
            collected.insert(collected.end(), generated.begin(), generated.end());
        }

        return collected;
    };

    // Replace duplicate phase names with unique names.
    auto fix_duplicate_phase_names = [](Vec<GeneralPhase>& phases)
    {
        Strings phasenames = vectorize(phases, RKT_LAMBDA(x, x.name()));
        phasenames = makeunique(phasenames, "!");
        auto i = 0;
        for(auto& phase : phases)
            phase.setName(phasenames[i++]);
    };

    Strings symbols = collect_all_element_symbols();

    Vec<GeneralPhase> allgeneralphases = collect_all_general_phases(symbols);

    fix_duplicate_phase_names(allgeneralphases);

    Vec<Phase> phases;
    phases.reserve(allgeneralphases.size());
    for(const auto& generalphase : allgeneralphases)
        phases.push_back(generalphase.convert(db, symbols));

    return phases;
}

Phases::operator Vec<Phase>() const
{
    return convert();
}

} // namespace Reaktoro
