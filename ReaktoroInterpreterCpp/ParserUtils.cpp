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

#include "ParserUtils.hpp"

// C++ includes
#include <istream>
#include <string>
#include <sstream>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <ReaktoroInterpreterCpp/Interpreter.hpp>

namespace Reaktoro {
namespace {

// Return a line with added `- ` right before the first non-space char in the line.
auto fixHyphen(std::string line) -> std::string
{
    // Get the index of first char that is not a space
    auto i = line.find_first_not_of(" ");

    // Check if the line is a comment line, an empty line, or a line with spaces only
    if(i == std::string::npos || line[i] == '#')
        return line;

    // Add `- ` right before the first non-space char in the line
    return line.substr(0, i) + "- " + line.substr(i);
}

/// Return a line with `Mixture` followed by `|` if needed.
auto fixMixture(std::string line) -> std::string
{
    auto i = line.find("Mixture:");

    if(i == std::string::npos)
        return line; // there is no `Mixture:` in the line

    auto j = line.find('#');

    if(j < i)
        return line; // `Mixture:` appears after comment symbol #

    auto k = line.find_first_not_of(" ", i + 8);

    if(k < j) // `Mixture:` is followed by an inline mixture recipe
        return line;

    // Return modified line with `Mixture: |`
    return line.substr(0, i + 8) + " |" + line.substr(i + 8);
}

} // namespace

MixtureCompound::MixtureCompound()
{}

MixtureCompound::MixtureCompound(std::string compound)
{
    auto words = split(compound);
    Assert(words.size() == 3, "Could not parse `" + compound + "` into a MixtureCompound.",
        "Expecting three words representing a value, its units, and a corresponding entity name, e.g., 1 kg H2O");
    value = tofloat(words[0]);
    units = words[1];
    entity = words[2];
    Assert(value >= 0.0, "Could not parse `" + compound + "` into a MixtureCompound.",
        "Expecting a non-negative amount for compound `" + entity + "`.");
}

auto operator>>(const Node& node, std::string& x) -> void
{
    x = str(valnode(node));
}

auto operator>>(const Node& node, ValueUnits& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() == 2, "Could not parse `" + node + "` into a ValueUnits instance.",
        "Expecting two words representing the value and its units, e.g., 300 kelvin, 50 moles");
    x.value = tofloat(words[0]);
    x.units = words[1];
}

auto operator>>(const Node& node, EntityValueUnits& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() == 3, "Could not parse `" + node + "` into an EntityValueUnits instance.",
        "Expecting three words representing an entity name, its value, and its units, e.g., Calcite 100 g");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.units = words[2];
}

auto operator>>(const Node& node, MixtureCompound& x) -> void
{
    x = MixtureCompound(str(node));
}

auto operator>>(const Node& node, MixtureNode& x) -> void
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
        "Expecting either an inline mixture of compounds "
        "(e.g., Mixture: 1 kg H2O; 1 mmol NaCl), and their amounts, "
        "or this list of compounds in subsequent indented lines (e.g., "
        "\nMixture:\n  - 1 kg H2O\n  - 1 mmol NaCl");
}

auto operator>>(const Node& node, EquilibriumConstraintNode::pH& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 0, "Could not parse the `pH` constraint.",
        "Expecting at least a value for the pH of the solution.");
    x.value = tofloat(words[0]);
    x.titrant1 = words.size() > 1 ? words[1] : "";
    x.titrant2 = words.size() > 2 ? words[2] : "";
}

auto operator>>(const Node& node, EquilibriumConstraintNode::SpeciesAmount& x) -> void
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

auto operator>>(const Node& node, EquilibriumConstraintNode::SpeciesActivity& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 1, "Could not parse the node `" + node + "`.",
        "Expecting at least a species name and a value for the species activity, e.g., `Activity: O2(g) 0.2`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.titrant1 = words.size() > 2 ? words[2] : "";
    x.titrant2 = words.size() > 3 ? words[3] : "";
}

auto operator>>(const Node& node, EquilibriumConstraintNode::PhaseAmount& x) -> void
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

auto operator>>(const Node& node, EquilibriumConstraintNode::PhaseVolume& x) -> void
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

auto operator>>(const Node& node, PlotNode& x) -> void
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

auto operator>>(const Node& node, EquilibriumNode& x) -> void
{
    ProcessFunction process_temperature = [&](const Node& child)
        { child >> x.temperature; };

    ProcessFunction process_pressure = [&](const Node& child)
        { child >> x.pressure; };

    ProcessFunction process_mixture = [&](const Node& child)
        { child >> x.mixture; };

    ProcessFunction process_pH = [&](const Node& child)
        { EquilibriumConstraintNode::pH c; child >> c; x.pH.push_back(c); };

    ProcessFunction process_species_amounts = [&](const Node& child)
        { EquilibriumConstraintNode::SpeciesAmount c; child >> c; x.species_amounts.push_back(c); };

    ProcessFunction process_species_activities = [&](const Node& child)
        { EquilibriumConstraintNode::SpeciesActivity c; child >> c; x.species_activities.push_back(c); };

    ProcessFunction process_phase_amounts = [&](const Node& child)
        { EquilibriumConstraintNode::PhaseAmount c; child >> c; x.phase_amounts.push_back(c); };

    ProcessFunction process_phase_volumes = [&](const Node& child)
        { EquilibriumConstraintNode::PhaseVolume c; child >> c; x.phase_volumes.push_back(c); };

    ProcessFunction process_inert_species = [&](const Node& child)
        { EntityValueUnits c; child >> c; x.inert_species.push_back(c); };

    ProcessFunction process_inert_phases = [&](const Node& child)
        { x.inert_phases = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"temperature"     , process_temperature},
        {"pressure"        , process_pressure},
        {"mixture"         , process_mixture},
        {"ph"              , process_pH},
        {"amount"          , process_species_amounts},
        {"activity"        , process_species_activities},
        {"speciesamount"   , process_species_amounts},
        {"speciesactivity" , process_species_activities},
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

auto operator>>(const Node& node, KineticsNode& x) -> void
{
    ProcessFunction process_initial_condition = [&](const Node& child)
        { child >> x.initial_condition; };

    ProcessFunction process_duration = [&](const Node& child)
        { child >> x.duration; };

    ProcessFunction process_plot = [&](const Node& child)
        { PlotNode p; child >> p; x.plots.push_back(p); };

    ProcessFunction process_kinetic_species = [&](const Node& child)
        { x.kinetic_species = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"initialcondition", process_initial_condition},
        {"initialstate"    , process_initial_condition},
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

auto operator>>(const Node& node, MineralReactionNode& x) -> void
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

    // Initialize the identifier of the chemical state
    x.mineral = identifier(node);

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }

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

auto preprocess(std::string script) -> std::string
{
    std::istringstream iss(script);
    return preprocess(iss);
}

auto preprocess(std::istream& stream) -> std::string
{
    std::stringstream ss;
    std::string line;
    while(std::getline(stream, line))
    {
        line = fixHyphen(line);
        ss << line << std::endl;
    }
    return ss.str();
}

} // namespace Reaktoro
