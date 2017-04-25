// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <unsupported/cpp-interpreter/Operators.hpp>

namespace Reaktoro {
namespace kwd {

auto operator>>(const Node& node, std::string& x) -> void
{
    x = str(valnode(node));
}

auto operator>>(const Node& node, ValueUnits& x) -> void
{
    x = ValueUnits(str(valnode(node)));
}

auto operator>>(const Node& node, EntityValueUnits& x) -> void
{
    x = EntityValueUnits(str(valnode(node)));
}

auto operator>>(const Node& node, ValueUnitsEntity& x) -> void
{
    x = ValueUnitsEntity(str(valnode(node)));
}

auto operator>>(const Node& node, Recipe& x) -> void
{
    Node val = valnode(node);

    if(val.IsScalar())
    {
        auto words = split(str(val), ";\n");
        for(auto word : words)
            x.push_back(word);
    }
    else if(val.IsSequence())
    {
        for(auto child : val)
            x.push_back(str(child));
    }
    else RuntimeError("Could not parse node `" + node + "`.",
        "Expecting either an inline list of compounds and their amounts "
        "(e.g., `Recipe: 1 kg H2O; 1 mmol NaCl`) "
        "or a multi-line list of compounds and their amounts (e.g., "
        "\nRecipe:\n  1 kg H2O\n  1 mmol NaCl");
}

auto operator>>(const Node& node, pH& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 0, "Could not parse the `pH` constraint.",
        "Expecting at least a value for the pH of the solution.");
    x.value = tofloat(words[0]);
    x.titrant1 = words.size() > 1 ? words[1] : "";
    x.titrant2 = words.size() > 2 ? words[2] : "";
}

auto operator>>(const Node& node, SpeciesAmount& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 2, "Could not parse the node `" + node + "`.",
        "Expecting at least a species name, an amount, and its units, e.g., `Amount: CO2(g) 1 mol`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.units = words[2];
    x.titrant1 = words.size() > 3 ? words[3] : "";
    x.titrant2 = words.size() > 4 ? words[4] : "";
}

auto operator>>(const Node& node, SpeciesActivity& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 1, "Could not parse the node `" + node + "`.",
        "Expecting at least a species name and a value for the species activity, e.g., `Activity: O2(g) 0.2`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.titrant1 = words.size() > 2 ? words[2] : "";
    x.titrant2 = words.size() > 3 ? words[3] : "";
}

auto operator>>(const Node& node, SpeciesFugacity& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 2, "Could not parse the node `" + node + "`.",
        "Expecting at least a species name, a value for the species fugacity, and its pressure units, e.g., `Fugacity: CO2(g) 40 Pa`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.units = words[2];
    x.titrant1 = words.size() > 3 ? words[3] : "";
    x.titrant2 = words.size() > 4 ? words[4] : "";
}

auto operator>>(const Node& node, PhaseAmount& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 2, "Could not parse the node `" + node + "`.",
        "Expecting at least a phase name, an amount, and its units, e.g., `PhaseAmount: Calcite 100 g`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.units = words[2];
    x.titrant1 = words.size() > 3 ? words[3] : "";
    x.titrant2 = words.size() > 4 ? words[4] : "";
}

auto operator>>(const Node& node, PhaseVolume& x) -> void
{
    Node rhs = valnode(node);
    auto words = split(str(valnode(node)));
    Assert(words.size() > 2, "Could not parse the node `" + node + "`.",
        "Expecting at least a phase name, a volume value, and its units, e.g., `PhaseVolume: Calcite 0.5 m3`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.units = words[2];
    x.titrant1 = words.size() > 3 ? words[3] : "";
    x.titrant2 = words.size() > 4 ? words[4] : "";
}

auto operator>>(const Node& node, Plot& x) -> void
{
    ProcessFunction process_name    = [&](const Node& child) { child >> x.name; };
    ProcessFunction process_x       = [&](const Node& child) { child >> x.x; };
    ProcessFunction process_y       = [&](const Node& child) { child >> x.y; };
    ProcessFunction process_xlabel  = [&](const Node& child) { child >> x.xlabel; };
    ProcessFunction process_ylabel  = [&](const Node& child) { child >> x.ylabel; };
    ProcessFunction process_ytitles = [&](const Node& child) { child >> x.ytitles; };
    ProcessFunction process_key     = [&](const Node& child) { child >> x.key; };

    std::map<std::string, ProcessFunction> fmap = {
        {"name"    , process_name},
        {"x"       , process_x},
        {"y"       , process_y},
        {"xlabel"  , process_xlabel},
        {"ylabel"  , process_ylabel},
        {"ytitles" , process_ytitles},
        {"key"     , process_key},
    };

    // Initialize the identifier of the chemical state
    x.name = identifier(node);

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }
}

auto operator>>(const Node& node, EquilibriumProblem& x) -> void
{
    ProcessFunction process_temperature = [&](const Node& child)
        { child >> x.temperature; };

    ProcessFunction process_pressure = [&](const Node& child)
        { child >> x.pressure; };

    ProcessFunction process_recipe = [&](const Node& child)
        { child >> x.recipe; };

    ProcessFunction process_pH = [&](const Node& child)
        { pH c; child >> c; x.ph.push_back(c); };

    ProcessFunction process_species_amounts = [&](const Node& child)
        { SpeciesAmount c; child >> c; x.species_amounts.push_back(c); };

    ProcessFunction process_species_activities = [&](const Node& child)
        { SpeciesActivity c; child >> c; x.species_activities.push_back(c); };

    ProcessFunction process_species_fugacities = [&](const Node& child)
        { SpeciesFugacity c; child >> c; x.species_fugacities.push_back(c); };

    ProcessFunction process_phase_amounts = [&](const Node& child)
        { PhaseAmount c; child >> c; x.phase_amounts.push_back(c); };

    ProcessFunction process_phase_volumes = [&](const Node& child)
        { PhaseVolume c; child >> c; x.phase_volumes.push_back(c); };

    ProcessFunction process_inert_species = [&](const Node& child)
        { EntityValueUnits c; child >> c; x.inert_species.push_back(c); };

    ProcessFunction process_inert_phases = [&](const Node& child)
        { x.inert_phases = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"temperature"     , process_temperature},
        {"pressure"        , process_pressure},
        {"recipe"         , process_recipe},
        {"ph"              , process_pH},
        {"amount"          , process_species_amounts},
        {"activity"        , process_species_activities},
        {"fugacity"        , process_species_fugacities},
        {"speciesamount"   , process_species_amounts},
        {"speciesactivity" , process_species_activities},
        {"speciesfugacity" , process_species_fugacities},
        {"phaseamount"     , process_phase_amounts},
        {"phasevolume"     , process_phase_volumes},
        {"inertspecies"    , process_inert_species},
        {"inertphases"     , process_inert_phases},
    };

    // Initialize the identifier of the chemical state
    x.stateid = identifier(node);

    // Prevent an empty identifier for the chemical state
    if(x.stateid.empty()) x.stateid = "State";

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }
}

auto operator>>(const Node& node, EquilibriumPath& x) -> void
{
    ProcessFunction process_initial_state = [&](const Node& child)
        { child >> x.initial_state; };

    ProcessFunction process_final_state = [&](const Node& child)
        { child >> x.final_state; };

    ProcessFunction process_inert_species = [&](const Node& child)
        { x.inert_species = split(str(valnode(child)), " ;"); };

    ProcessFunction process_plot = [&](const Node& child)
        { Plot p; child >> p; x.plots.push_back(p); };

    std::map<std::string, ProcessFunction> fmap = {
        {"initialstate"    , process_initial_state},
        {"finalstate"      , process_final_state},
        {"inertspecies"    , process_inert_species},
        {"plot"            , process_plot},
    };

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }

    // Assert an initial state was provided.
    Assert(x.initial_state.size(), "Could not parse node `" + node + "`",
        "Expecting an `InitialState` keyword-value pair.")

    // Assert a final state was provided.
    Assert(x.final_state.size(), "Could not parse node `" + node + "`",
        "Expecting an `FinalState` keyword-value pair.")
}

auto operator>>(const Node& node, KineticPath& x) -> void
{
    ProcessFunction process_initial_condition = [&](const Node& child)
        { child >> x.initial_condition; };

    ProcessFunction process_duration = [&](const Node& child)
        { child >> x.duration; };

    ProcessFunction process_plot = [&](const Node& child)
        { Plot p; child >> p; x.plots.push_back(p); };

    ProcessFunction process_inert_species = [&](const Node& child)
        { x.inert_species = split(str(valnode(child)), " ;"); };

    ProcessFunction process_kinetic_species = [&](const Node& child)
        { x.kinetic_species = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"initialcondition", process_initial_condition},
        {"initialstate"    , process_initial_condition},
        {"inertspecies"    , process_inert_species},
        {"kineticspecies"  , process_kinetic_species},
        {"duration"        , process_duration},
        {"plot"            , process_plot},
    };

    // Initialize the identifier of the chemical state
    x.stateid = identifier(node);

    // Prevent an empty identifier for the chemical state
    if(x.stateid.empty()) x.stateid = "State";

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }

    // Assert an initial condition was provided.
    Assert(x.initial_condition.size(), "Could not parse node `" + node + "`",
        "Expecting an `InitialCondition` keyword-value pair.")
}

auto operator>>(const Node& node, MineralReaction& x) -> void
{
    ProcessFunction process_mineral = [&](const Node& child)
        { child >> x.mineral; };

    ProcessFunction process_equation = [&](const Node& child)
        { x.equation = str(valnode(child)); };

    ProcessFunction process_mechanism = [&](const Node& child)
        { x.mechanisms.push_back(str(valnode(child))); };

    ProcessFunction process_ssa = [&](const Node& child)
        { child >> x.ssa; };

    std::map<std::string, ProcessFunction> fmap = {
        {"mineral"             , process_mineral},
        {"equation"            , process_equation},
        {"mechanism"           , process_mechanism},
        {"specificsurfacearea" , process_ssa},
        {"ssa"                 , process_ssa},
    };

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }

    // Initialize the identifier of the mineral reaction if still empty
    x.mineral = x.mineral.empty() ? identifier(node) : x.mineral;

    // Assert a mineral name was given
    Assert(x.mineral.size(), "Could not parse node `" + node + "`",
        "Expecting a `Mineral` keyword-value pair (e.g., `Mineral: Quartz`.")

    // Assert a reaction equation was given
    Assert(x.mineral.size(), "Could not parse node `" + node + "`",
        "Expecting a `Equation` keyword-value pair (e.g., `Equation: Calcite = Ca++ + CO3--`.")

    // Assert a reaction equation was given
    Assert(x.ssa.value > 0, "Could not parse node `" + node + "`",
        "Expecting a `SpecificSurfaceArea` keyword-value pair (e.g., `SpecificSurfaceArea: 10 cm2/g`.")

    // Assert at least one mineral kinetic mechanism was given
    Assert(x.mechanisms.size(), "Could not parse node `" + node + "`",
        "Expecting at least one `Mechanism` keyword-value pair (e.g., `Mechanism: logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0`.")
}

auto operator>>(const Node& node, Concentrations& x) -> void
{
    Node val = valnode(node);

    if(val.IsScalar())
    {
        auto words = split(str(val), ";\n");
        for(auto word : words)
            x.push_back(word);
    }
    else if(val.IsSequence())
    {
        for(auto child : val)
            x.push_back(str(child));
    }
    else RuntimeError("Could not parse node `" + node + "`.",
        "Expecting either an inline list of compounds and their concentrations "
        "(e.g., `Concentrations: Na 203; Ca 87`) "
        "or a multi-line list of compounds and their concentrations (e.g., "
        "\nConcentrations:\n  Na 203\n  Ca  87");
}

auto operator>>(const Node& node, SpeciationProblem& x) -> void
{
    ProcessFunction process_temperature = [&](const Node& child)
        { child >> x.temperature; };

    ProcessFunction process_pressure = [&](const Node& child)
        { child >> x.pressure; };

    ProcessFunction process_concentrations = [&](const Node& child)
        { child >> x.concentrations; };

    ProcessFunction process_pH = [&](const Node& child)
        { pH c; child >> c; x.ph.push_back(c); };

    ProcessFunction process_species_activities = [&](const Node& child)
        { SpeciesActivity c; child >> c; x.species_activities.push_back(c); };

    ProcessFunction process_species_fugacities = [&](const Node& child)
        { SpeciesFugacity c; child >> c; x.species_fugacities.push_back(c); };

    ProcessFunction process_inert_species = [&](const Node& child)
        { EntityValueUnits c; child >> c; x.inert_species.push_back(c); };

    ProcessFunction process_inert_phases = [&](const Node& child)
        { x.inert_phases = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"temperature"     , process_temperature},
        {"pressure"        , process_pressure},
        {"concentrations"  , process_concentrations},
        {"ph"              , process_pH},
        {"activity"        , process_species_activities},
        {"fugacity"        , process_species_fugacities},
        {"speciesactivity" , process_species_activities},
        {"speciesfugacity" , process_species_fugacities},
        {"inertspecies"    , process_inert_species},
        {"inertphases"     , process_inert_phases},
    };

    // Initialize the identifier of the chemical state
    x.stateid = identifier(node);

    // Prevent an empty identifier for the chemical state
    if(x.stateid.empty()) x.stateid = "Speciation";

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }
}

auto operator>>(const Node& node, PhreeqcKeyword& x) -> void
{
    ProcessFunction process_database = [&](const Node& child)
	{
    	child >> x.database;
	};

    ProcessFunction process_input = [&](const Node& child)
	{
		child >> x.input;
    };

    ProcessFunction process_output = [&](const Node& child)
	{
		child >> x.output;
    };

    std::map<std::string, ProcessFunction> fmap = {
        {"database" , process_database},
        {"input"    , process_input},
        {"output"   , process_output},
    };

    // The right-hand side node of `node`
    Node val = valnode(node);

    // Initialize the identifier of the chemical state
    x.stateid = identifier(node);

    // Prevent an empty identifier for the chemical state
    if(x.stateid.empty()) x.stateid = "StatePhreeqc";

	// Otherwise, process each child of the PHREEQC keyword
	for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }

	// Assert both database and input has been given
	Assert(x.database.size(), "Could not parse node `" + node + "`.",
		"Expecting a database file name, e.g., `Database: phreeqc.dat`");

	Assert(x.input.size(), "Could not parse node `" + node + "`.",
		"Expecting an input script, e.g., `Input: input.dat`");
}

} // namespace kwd
} // namespace Reaktoro
