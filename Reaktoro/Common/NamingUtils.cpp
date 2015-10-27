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

#include "NamingUtils.hpp"

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>

namespace Reaktoro {
namespace {

/// The alternative names for water
std::vector<std::string> alternative_water_names = {"H2O(l)", "H2O", "H2O@", "H2O(aq)"};

/// The map with alternative names for charged species
std::map<std::pair<std::string, double>, std::vector<std::string>> alternative_charged_names;

/// The map with alternative names for neutral species
std::map<std::string, std::vector<std::string>> alternative_neutral_names;

} // namespace

auto alternativeWaterNames() -> std::vector<std::string>&
{
    return alternative_water_names;
}

auto alternativeChargedSpeciesNames(std::string name) -> std::vector<std::string>&
{
    // Get the charge out of the species name (Ca++ => +2, CO3-- => -2)
    const int charge = chargeInSpeciesName(name);

    // Check if the species has charge
    Assert(charge != 0, "Could not get an alternative charged species name.",
        "The given species `" + name + "` is not charged or its name does not "
        "follow the convention `H+`, `Ca++`, `HCO3-`, `SO4--`, and so forth.");

    // Get the base name of the charged species
    const std::string base = baseNameChargedSpecies(name);

    // Auxiliary variables
    const char sign = (charge < 0) ? '-' : '+';
    const int abs_charge = std::abs(charge);
    const std::string str_charge = sign + (abs_charge == 1 ? "" : std::to_string(abs_charge));

    // Create alternative names for the given charged species
    std::vector<std::string> alternatives;
    alternatives.push_back(base + std::string(abs_charge, sign)); // e.g.: Ca++
    alternatives.push_back(base + str_charge);                    // e.g.: Ca+2
    alternatives.push_back(base + "[" + str_charge + "]");        // e.g.: Ca[+2]

    auto res = alternative_charged_names.insert({{base, charge}, unique(alternatives)});

    return res.first->second; // return the iterator to the existing or just added alternative names
}

auto alternativeNeutralSpeciesNames(std::string name) -> std::vector<std::string>&
{
    // Get the base name of the neutral species
    const std::string base = baseNameNeutralSpecies(name);

    std::vector<std::string> alternatives;
    alternatives.push_back(base + "(aq)"); // e.g.: CO2(aq)
    alternatives.push_back(base);          // e.g.: CO2
    alternatives.push_back(base + "@");    // e.g.: CO2@
    alternatives.push_back(base + ",aq");  // e.g.: CO2,aq

    auto res = alternative_neutral_names.insert({base, unique(alternatives)});

    return res.first->second; // return the iterator to the existing or just added alternative names
}

auto conventionalWaterName() -> std::string
{
    return "H2O(l)";
}

auto conventionalChargedSpeciesName(std::string name) -> std::string
{
    // Get the last char of the string
    const char last = *name.rbegin();

    // Check for Ca[2+], H[+], CO3[2-], Cl[-] convention
    if(last == )
    // Check for Ca++, H+, CO3--, Cl- convention
    if(name.rbegin())
    auto alternativeChargedSpeciesNames(name);
}

auto conventionalNeutralSpeciesName(std::string name) -> std::string
{

}

auto isAlternativeWaterName(std::string trial) -> bool
{
    return contained(trial, alternativeWaterNames());
}

auto isAlternativeChargedSpeciesName(std::string trial, std::string name) -> bool
{
    return contained(trial, alternativeChargedSpeciesNames(name));
}

auto isAlternativeNeutralSpeciesName(std::string trial, std::string name) -> bool
{
    return contained(trial, alternativeNeutralSpeciesNames(name));
}

auto baseNameChargedSpecies(std::string name) -> std::string
{
    // Get the sign of the charge by inspecting last char in the string
    const char sign = *name.rbegin();

    // Check if the species name is neutral
    if(sign != '-' && sign != '+') return name;

    // Find the position of the first char that matches sign
    const auto pos = name.find(sign);

    return name.substr(0, pos);
}

auto baseNameNeutralSpecies(std::string name) -> std::string
{
    // Find the position of the first char that matches sign
    const auto pos = name.find("(aq)");

    // Check if the species name has the suffix `(aq)`
    if(pos == std::string::npos)
        return name;

    return name.substr(0, pos);
}

auto chargeInSpeciesName(std::string name) -> double
{
    // Get the sign of the charge by inspecting last char in the string
    const char sign = *name.rbegin();

    // Check if the species name is neutral
    if(sign != '-' && sign != '+') return 0;

    // Find the position of the first char that matches sign
    const auto pos1 = name.find(sign);

    // Find the position of the first char that does not match sign after `pos1`
    const auto pos2 = name.find_first_not_of(sign, pos1);

    // Check for invalid cases such as `Ca-+`
    Assert(pos2 == std::string::npos,
        "Could not read the charge of the species from its name.",
        "The given species name `" + name + "` does not have a "
        "continuous trailing sequence of signs `-` or `+`.");

    const double factor = (sign == '-') ? -1.0 : +1.0;

    return factor * (name.size() - pos1);
}

} // namespace Reaktoro
