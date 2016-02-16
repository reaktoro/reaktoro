// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Processors.hpp"

// Reaktoro includes
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>

// Interpreter includes
#include "Initializers.hpp"
#include "Interpreter.hpp"
#include "Keywords.hpp"
#include "Operators.hpp"
#include "Yaml.hpp"

namespace Reaktoro {
namespace {

/// Return all entity names it can find in list of EntityValueUnits objects.
template<typename Triplet>
auto collectEntities(std::vector<std::string>& list, const std::vector<Triplet>& triplets) -> void
{
    for(const Triplet& t : triplets)
        list.push_back(t.entity);
}

/// Return all titrant names it can find in a list of EquilibriumConstraintNode objects.
template<typename Constraint>
auto collectTitrants(std::vector<std::string>& list, const std::vector<Constraint>& constraints) -> void
{
    for(const Constraint& c : constraints)
    {
        list.push_back(c.entity);
        list.push_back(c.titrant1);
        list.push_back(c.titrant2);
    }
}

/// Return all compounds (e.g., species and titrant names) it can find in an Equilibrium instance.
/// This method is used to setup an automatic chemical system composed of only aqueous species.
auto collectCompounds(const kwd::EquilibriumProblem& e) -> std::vector<std::string>
{
    std::vector<std::string> list;
    collectEntities(list, e.recipe);
    collectTitrants(list, e.ph);
    collectTitrants(list, e.species_amounts);
    collectTitrants(list, e.species_activities);
    collectTitrants(list, e.species_fugacities);
    collectEntities(list, e.inert_species);
    remove(list, [](std::string x) { return x.empty(); });
    return list;
}

/// Return a list of compounds names that are present as mineral species in a database.
auto filterMineralSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::vector<std::string> minerals;
    for(auto compound : compounds)
        if(database.containsMineralSpecies(compound))
            minerals.push_back(compound);
    return minerals;
}

/// Return a list of compounds names that are present as gaseous species in a database.
auto filterGaseousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::vector<std::string> gases;
    for(auto compound : compounds)
        if(database.containsGaseousSpecies(compound))
            gases.push_back(compound);
    return gases;
}

} // namespace

auto processMineralReactionNode(InterpreterState& istate, const Node& node) -> void
{
    // Convert the yaml node into a keyword
    kwd::MineralReaction keyword; node >> keyword;

    // Convert the keyword into a MineralReaction instance
    MineralReaction reaction;
    initializeMineralReaction(reaction, keyword);

    // Add the new mineral reaction to the list
    istate.mineral_reactions.push_back(reaction);
}

auto processEquilibriumNode(InterpreterState& istate, const Node& node) -> void
{
    // Convert the yaml node into a keyword
    kwd::EquilibriumProblem keyword; node >> keyword;

    // Check if a chemical system needs to be defined
    if(istate.system.numSpecies() == 0)
    {
        // Collect all compounds that appears in the definition of the equilibrium problem
        auto compounds = collectCompounds(keyword);

        // Determine if there are gaseous and mineral species among those compounds
        auto mineral_species = filterMineralSpecies(compounds, istate.database);
        auto gaseous_species = filterGaseousSpecies(compounds, istate.database);

        // Add the aqueous phase using the compounds as the initializer
        istate.editor.addAqueousPhaseWithCompounds(compounds);

        // Add a gaseous phase if there were gaseous species among the compound names
        if(gaseous_species.size())
            istate.editor.addGaseousPhaseWithSpecies(gaseous_species);

        // Add a mineral phase for each mineral species among the compound names
        for(auto x : mineral_species)
            istate.editor.addMineralPhaseWithSpecies({x});

        // Initialize the chemical system
        istate.system = istate.editor;
    }

    // Initialize the equilibrium problem using the just initialized chemical system
    EquilibriumProblem problem(istate.system);
    initializeEquilibriumProblem(problem, keyword);

    // Initialize the chemical state
    ChemicalState state(istate.system);

    // Initialize the amounts of the inert species
    for(auto s : keyword.inert_species)
        state.setSpeciesAmount(s.entity, s.value, s.units);

    // Perform the equilibrium calculation
    equilibrate(state, problem);

    // Output the resulting chemical state
    state.output(keyword.stateid + ".dat");

    // Store the resulting chemical state with given state ID
    istate.states[keyword.stateid] = state;
}

auto processEquilibriumPathNode(InterpreterState& istate, const Node& node) -> void
{
    // Convert the yaml node into a kinetic keyword
    kwd::EquilibriumPath keyword; node >> keyword;

    // Initialize the kinetic path instance
    EquilibriumPath path(istate.system);
    initializeEquilibriumPath(path, keyword);

    // Alias to the initial and final states
    auto& statei = istate.states[keyword.initial_state];
    auto& statef = istate.states[keyword.final_state];

    // Solve the equilibrium path problem
    path.solve(statei, statef);
}

auto processKineticPathNode(InterpreterState& istate, const Node& node) -> void
{
    // Convert the yaml node into a kinetic keyword
    kwd::KineticPath keyword; node >> keyword;

    // Initialize the ReactionSystem instance
    for(auto r : istate.mineral_reactions)
        istate.editor.addMineralReaction(r);

    // Set the reactions of the chemical system
    istate.reactions = istate.editor;

    // Initialize the kinetic path instance
    KineticPath path(istate.reactions);
    initializeKineticPath(path, keyword);

    // Alias to the initial condition state
    auto& state = istate.states[keyword.initial_condition];

    // The duration time and its units
    const auto duration = keyword.duration.value;
    const auto units = keyword.duration.units;

    // Solve the kinetics problem
    path.solve(state, 0, duration, units);

    // Output the final chemical state with given state id
    state.output(keyword.stateid + ".dat");

    // Store the final chemical state into the map of states
    istate.states[keyword.stateid] = state;
}

} // namespace Reaktoro
