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

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>
#include <Reaktoro/Thermodynamics/Phases/AqueousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// Reaktoro includes
#include "Initializers.hpp"
#include "Interpreter.hpp"
#include "Keywords.hpp"
#include "Operators.hpp"
#include "Utils.hpp"
#include "Yaml.hpp"

namespace Reaktoro {

auto processDatabaseNode(InterpreterState& istate, const Node& node) -> void
{
    istate.database = Database(str(valnode(node)));
    istate.editor = ChemicalEditor(istate.database);
    istate.elements = identifyElements(istate.compounds, istate.database);
}

auto processAqueousPhaseNode(InterpreterState& istate, const Node& node) -> void
{
    const std::string phasename = identifier(node);
    Assert(phasename.empty(), "Could not set the aqueous phase with `" + node + "`.",
        "The name of the aqueous phase is `Aqueous` and cannot be changed.");
    auto compounds = str(valnode(node));
    if(compounds.find("auto") != std::string::npos)
        istate.editor.addAqueousPhaseWithElements(istate.elements);
    else istate.editor.addAqueousPhase(compounds);
}

auto processGaseousPhaseNode(InterpreterState& istate, const Node& node) -> void
{
    const std::string phasename = identifier(node);
    Assert(phasename.empty(), "Could not set the gaseous phase with `" + node + "`.",
        "The name of the gaseous phase is `Gaseous` and cannot be changed.");
    auto compounds = str(valnode(node));
    if(compounds.find("auto") != std::string::npos)
        istate.editor.addGaseousPhaseWithElements(istate.elements);
    else istate.editor.addGaseousPhase(compounds);
}

auto processMineralPhaseNode(InterpreterState& istate, const Node& node) -> void
{
    const std::string phasename = identifier(node);
    auto& phase = istate.editor.addMineralPhase(str(valnode(node)));
    if(!phasename.empty())
        phase.setName(phasename);
}

auto processMineralsNode(InterpreterState& istate, const Node& node) -> void
{
    auto minerals = str(valnode(node));
    if(minerals.find("auto") != std::string::npos)
        for(auto mineral : istate.database.mineralSpeciesWithElements(istate.elements))
            istate.editor.addPhase(MineralPhase(mineral));
    else for(auto mineral : split(minerals, " ;"))
        istate.editor.addMineralPhase(mineral);
}

auto processChemicalModelNode(InterpreterState& istate, const Node& node) -> void
{
}

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

    // Initialize the equilibrium problem using the just initialized chemical system
    EquilibriumProblem problem(istate.system);
    initializeEquilibriumProblem(problem, keyword);

    // Initialize the chemical state
    EquilibriumState state(istate.system);

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

auto processPhreeqcNode(InterpreterState& istate, const Node& node) -> void
{
    // Convert the yaml node into a kinetic keyword
    kwd::PhreeqcKeyword keyword; node >> keyword;

    // Initialize a Phreeqc instance with given database and input script
    Phreeqc phreeqc;
    phreeqc.load(keyword.database);
    phreeqc.execute(keyword.input, keyword.output);

    // Create a chemical state instance
    ChemicalState state = phreeqc;

    // Initialize the chemical system
    istate.system = state.system();

    // Output the final chemical state with given state id
    state.output(keyword.stateid + ".dat");

    // Store the final chemical state into the map of states
    istate.states[keyword.stateid] = state;
}

} // namespace Reaktoro
