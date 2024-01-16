// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "NamingUtils.hpp"

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

#include <cctype>

namespace Reaktoro {
namespace {

/// The alternative names for water
std::vector<std::string> alternative_water_names = {"H2O(l)", "H2O", "H2O@", "H2O(aq)"};

/// The map with alternative names for charged species
std::map<std::pair<std::string, double>, std::vector<std::string>> alternative_charged_names;

/// The map with alternative names for neutral species
std::map<std::string, std::vector<std::string>> alternative_neutral_names;

} // namespace

auto alternativeWaterNames() -> std::vector<std::string>
{
    return alternative_water_names;
}

auto alternativeChargedSpeciesNames(std::string name) -> std::vector<std::string>&
{
    // Get the base name of the charged species and its charge
    const auto pair = splitChargedSpeciesName(name);
    const std::string base = pair.first;
    const int charge = pair.second;
    const int abs_charge = std::abs(charge);
    const char sign = (charge < 0) ? '-' : '+';
    const std::string sign_charge = sign + (abs_charge == 1 ? "" : std::to_string(abs_charge));
    const std::string charge_sign = (abs_charge == 1 ? "" : std::to_string(abs_charge)) + sign;

    // Check if the species has charge
    Assert(charge != 0, "Could not get an alternative charged species name.",
        "The given species `" + name + "` is not charged or its name does not "
        "follow the convention `H+`, `Ca++`, `HCO3-`, `SO4--`, and so forth.");

    // Create alternative names for the given charged species
    std::vector<std::string> alternatives;
    alternatives.push_back(base + std::string(abs_charge, sign)); // e.g.: Ca++
    alternatives.push_back(base + sign_charge);                   // e.g.: Ca+2
    alternatives.push_back(base + "[" + charge_sign + "]");       // e.g.: Ca[2+]

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
    auto pair = splitChargedSpeciesName(name);
    const std::string base = pair.first;
    const int charge = pair.second;
    const char sign = (charge < 0) ? '-' : '+';
    return base + std::string(std::abs(charge), sign);
}

auto conventionalNeutralSpeciesName(std::string name) -> std::string
{
    const std::string base = baseNameNeutralSpecies(name);
    return base + "(aq)";
}

auto isAlternativeWaterName(std::string trial) -> bool
{
    return contains(alternativeWaterNames(), trial);
}

auto isAlternativeChargedSpeciesName(std::string trial, std::string name) -> bool
{
    return contains(alternativeChargedSpeciesNames(name), trial);
}

auto isAlternativeNeutralSpeciesName(std::string trial, std::string name) -> bool
{
    return contains(alternativeNeutralSpeciesNames(name), trial);
}

auto splitChargedSpeciesName(std::string name) -> std::pair<std::string, double>
{
    const auto mid = name.find_first_of("-+[");

    const std::string base = name.substr(0, mid);

    // Check if the charged species name is of the form Ca[2+], H[+], Cl[-]
    if(name[mid] == '[')
    {
        const auto pos1 = mid;
        const auto pos2 = name.find(']', pos1);
        Assert(pos2 != std::string::npos, "Could not parse the species name `" + name + "`.",
            "The species name is missing the closing bracket `]`");
        const std::string inside = name.substr(pos1 + 1, pos2 - pos1 - 1);
        const char sign = name[pos2 - 1];
        const std::string numb = inside.size() > 1 ? inside.substr(0, inside.size() - 1) : "1";
        const int charge = (sign == '-') ? -std::stoi(numb) : std::stoi(numb);
        return {base, charge};
    }
    // Check if the charged species name is of the form Cl-, CO3--, SO4--, Ca++, Na+, H+
    if((name[mid] == '-' || name[mid] == '+') && name.back() == name[mid])
    {
        const int abscharge = name.size() - mid;
        const int charge = name[mid] == '-' ? -abscharge : +abscharge;
        return {base, charge};
    }
    // Check if the charged species name is of the form Cl-, CO3--, SO4--
    if((name[mid] == '-' || name[mid] == '+') && std::isdigit(name.back()))
    {
        const int charge = std::stoi(name.substr(mid));
        return {base, charge};
    }

    RuntimeError("Could not parse the species name `" + name + "`.",
        "The species name has no recognized format such as Ca++, Ca+2, or Ca[2+].");
}

auto baseNameChargedSpecies(std::string name) -> std::string
{
    return splitChargedSpeciesName(name).first;
}

auto baseNameNeutralSpecies(std::string name) -> std::string
{
    // Check the convention CO2(aq), CaCO3(aq)
    auto pos = name.find("(aq)");
    if(pos < name.size()) return name.substr(0, pos);

    // Check the convention CO2@, CaCO3@
    pos = name.find("@");
    if(pos < name.size()) return name.substr(0, pos);

    // Check the convention CO2,aq, CaCO3,aq
    pos = name.find(",aq");
    if(pos < name.size()) return name.substr(0, pos);

    // Cover the case where no suffix is used or the used suffix is unrecognized
    return name;
}

auto chargeInSpeciesName(std::string name) -> double
{
    return splitChargedSpeciesName(name).second;
}

auto splitSpeciesNameSuffix(std::string name) -> std::pair<std::string, std::string>
{
    if(name.back() != ')')
        return { name, "" };

    const auto pos = name.rfind('(');
    const auto suffix = name.substr(pos + 1, name.size() - pos - 2); // remove both ( and )

    for(auto i = 0; i < suffix.size(); ++i)
        if(isupper(suffix[i])) // if there is upper case char, then it is not species name suffix
            return { name, "" }; // it could be part of the chemical formula, e.g., Fe(CO3)

    return { name.substr(0, pos), suffix };
}

} // namespace Reaktoro
