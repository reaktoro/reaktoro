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

/// The type used to represent a node processor function
using ProcessFunction = std::function<void(const Node&)>;

/// Return a joined string with a node string representation.
auto operator+(std::string str, const Node& node) -> std::string
{
    std::stringstream res;
    res << str << node;
    return res.str();
}

/// Return a joined string with a node string representation.
auto operator+(const Node& node, std::string str) -> std::string
{
    std::stringstream res;
    res << node << str;
    return res.str();
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

auto operator>>(const Node& node, Mixture& x) -> void
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
    else RuntimeError("Could not parse node `" + node + "` into a Mixture.",
        "Expecting either an inline mixture of compounds "
        "(e.g., Mixture: 1 kg H2O; 1 mmol NaCl), and their amounts, "
        "or this list of compounds in subsequent indented lines (e.g., "
        "\nMixture:\n  - 1 kg H2O\n  - 1 mmol NaCl");
}

auto operator>>(const Node& node, EquilibriumConstraint::pH& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 0, "Could not parse the `pH` constraint.",
        "Expecting at least a value for the pH of the solution.");
    x.value = tofloat(words[0]);
    x.titrant1 = words.size() > 1 ? words[1] : "";
    x.titrant2 = words.size() > 2 ? words[2] : "";
}

auto operator>>(const Node& node, EquilibriumConstraint::SpeciesAmount& x) -> void
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

auto operator>>(const Node& node, EquilibriumConstraint::SpeciesActivity& x) -> void
{
    auto words = split(str(valnode(node)));
    Assert(words.size() > 1, "Could not parse the node `" + node + "`.",
        "Expecting at least a species name and a value for the species activity, e.g., `Activity: O2(g) 0.2`.");
    x.entity = words[0];
    x.value = tofloat(words[1]);
    x.titrant1 = words.size() > 2 ? words[2] : "";
    x.titrant2 = words.size() > 3 ? words[3] : "";
}

auto operator>>(const Node& node, EquilibriumConstraint::PhaseAmount& x) -> void
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

auto operator>>(const Node& node, EquilibriumConstraint::PhaseVolume& x) -> void
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

auto operator>>(const Node& node, Equilibrium& x) -> void
{
    ProcessFunction process_temperature = [&](const Node& child)
        { child >> x.temperature; };
        
    ProcessFunction process_pressure = [&](const Node& child)
        { child >> x.pressure; };
        
    ProcessFunction process_mixture = [&](const Node& child)
        { child >> x.mixture; };

    ProcessFunction process_pH = [&](const Node& child)
        { EquilibriumConstraint::pH c; child >> c; x.pH.push_back(c); };
        
    ProcessFunction process_species_amounts = [&](const Node& child)
        { EquilibriumConstraint::SpeciesAmount c; child >> c; x.species_amounts.push_back(c); };
        
    ProcessFunction process_species_activities = [&](const Node& child)
        { EquilibriumConstraint::SpeciesActivity c; child >> c; x.species_activities.push_back(c); };
        
    ProcessFunction process_phase_amounts = [&](const Node& child)
        { EquilibriumConstraint::PhaseAmount c; child >> c; x.phase_amounts.push_back(c); };
        
    ProcessFunction process_phase_volumes = [&](const Node& child)
        { EquilibriumConstraint::PhaseVolume c; child >> c; x.phase_volumes.push_back(c); };

    ProcessFunction process_inert_species = [&](const Node& child)
        { EntityValueUnits c; child >> c; x.inert_species.push_back(c); };

    ProcessFunction process_inert_phases = [&](const Node& child)
        { x.inert_phases = split(str(valnode(child)), " ;"); };

    std::map<std::string, ProcessFunction> fmap = {
        {"Temperature"     , process_temperature},
        {"Pressure"        , process_pressure},
        {"Mixture"         , process_mixture},
        {"pH"         ,      process_pH},
        {"Amount"          , process_species_amounts},
        {"Activity"        , process_species_activities},
        {"SpeciesAmount"   , process_species_amounts},
        {"SpeciesActivity" , process_species_activities},
        {"PhaseAmount"     , process_phase_amounts},
        {"PhaseVolume"     , process_phase_volumes},
        {"InertSpecies"    , process_inert_species},
        {"InertPhases"     , process_inert_phases},
    };

    // Initialize the identifier of the chemical state
    x.stateid = identifier(node);

    // Prevent an empty identifier for the chemical state
    if(x.stateid.empty()) x.stateid = "State";

    Node val = valnode(node);

    for(auto child : val)
    {
        std::string key = keyword(child);
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }
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

auto str(const Node& node) -> std::string
{
    std::stringstream ss; ss << node;
    return ss.str();
}

auto keynode(const Node& node) -> Node
{
    Assert(node.IsMap(), "Could not get the key of node `" + str(node) + "`.",
        "This is only allowed for nodes that are maps.");
    return node.begin()->first;
}

auto valnode(const Node& node) -> Node
{
    Assert(node.IsMap(), "Could not get the value of node `" + str(node) + "`.",
        "This is only allowed for nodes that are maps.");
    return node.begin()->second;
}

auto keyword(const Node& node) -> std::string
{
    std::string key = str(keynode(node));
    auto words = split(key);
    return words[0];
}

auto identifier(const Node& node) -> std::string
{
    std::string key = str(keynode(node));
    auto words = split(key);
    return words.size() > 1 ? words[1] : "";
}

} // namespace Reaktoro
