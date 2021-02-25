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

#pragma once

// C++ includes
#include <string>
#include <vector>
#include <cctype>

namespace Reaktoro {

/// Return a collection of alternative names to the given species name.
/// @param name The name of the species.
auto alternativeWaterNames() -> std::vector<std::string>;

/// Return a collection of alternative names to the given charged species name.
/// The argument `name` must follow the naming convention `H+`, `Ca++`, `Cl-`,
/// `HCO3-`, `CO3--`, and so forth. That is it, the charge is given as a suffix
/// with as much `-` or `+` as the number of charges.
/// @param name The name of the charged species.
auto alternativeChargedSpeciesNames(std::string name) -> std::vector<std::string>&;

/// Return a collection of alternative names to the given neutral species name.
/// The argument `name` must follow the naming convention `CO2(aq)`, `CaCO3(aq)`,
/// `NaCl(aq)`, and so forth. That is, the neutral species name has the suffix `(aq)`.
/// @param name The name of the neutral species.
auto alternativeNeutralSpeciesNames(std::string name) -> std::vector<std::string>&;

/// Return the conventional water name, which is `H2O(l)`.
auto conventionalWaterName() -> std::string;

/// Return the conventional charged species name adopted in Reaktoro.
/// This method will return the conventional charged species name
/// adopted in Reaktoro. For example, `Ca+2` results in `Ca++`, `CO3-2`
/// results in `CO3--`, `Mg[+2]` results in `Mg++`, and so forth.
/// @param name The name of a charged species.
auto conventionalChargedSpeciesName(std::string name) -> std::string;

/// Return the conventional neutral species name adopted in Reaktoro.
/// This method will return the conventional neutral species name
/// adopted in Reaktoro. For example, `CO2` results in `CO2(aq)`, `CaCO3@`
/// results in `CaCO3(aq)`, `NaCl,aq` results in `NaCl(aq)`, and so forth.
/// @param name The name of a charged species.
auto conventionalNeutralSpeciesName(std::string name) -> std::string;

/// Return true if a `trial` name is an alternative to a water species name.
/// @param trial The trial name that is being checked as an alternative to `H2O(l)`
auto isAlternativeWaterName(std::string trial) -> bool;

/// Return true if a `trial` name is an alternative to a given charged species `name`.
/// @param trial The trial name that is being checked as an alternative to `name`
/// @param name The name of the charged species with convention `H+`, `Ca++`, `CO3--`, etc.
auto isAlternativeChargedSpeciesName(std::string trial, std::string name) -> bool;

/// Return true if a `trial` name is an alternative to a given neutral species `name`.
/// @param trial The trial name that is being checked as an alternative to `name`
/// @param name The name of the neutral species with convention `CO2(aq)`, `CaCO3(aq)`.
auto isAlternativeNeutralSpeciesName(std::string trial, std::string name) -> bool;

/// Return a pair with the base name of a charged species and its electrical charge.
/// The name of the charge species has to have the following naming conventions:
/// Ca++, Ca+2, and Ca[2+], or Na+, Na[+], or CO3--, CO3-2, and CO3[2-].
/// @param name The name of the charged species.
auto splitChargedSpeciesName(std::string name) -> std::pair<std::string, double>;

/// Return the name of a charged species without charge suffix.
/// The name of the charge species has to have the following naming conventions:
/// Ca++, Ca+2, and Ca[2+], or Na+, Na[+], or CO3--, CO3-2, and CO3[2-].
/// @param name The name of the charged species.
auto baseNameChargedSpecies(std::string name) -> std::string;

/// Return the name of a neutral species without suffix `(aq)`.
/// @param name The name of the charged species.
auto baseNameNeutralSpecies(std::string name) -> std::string;

/// Return the electrical charge in a species name.
/// The charge in the species name must be represended as
/// `H+`, `Ca++`, `HCO3-`, `CO3--`, `SO4--`, and so forth.
/// @param name The name of the species
auto chargeInSpeciesName(std::string name) -> double;

/// Split name and suffix from a substance name or chemical formula.
/// The suffix must start with `(`, end with `)` and cannot contain upper case
/// characters.
///
/// Example:
/// ~~~
/// using namespace Reaktoro;
/// const auto [name, suffix] = splitSpeciesNameSuffix("H2O(aq)");           // name: "H2O",     suffix: "aq"
/// const auto [name, suffix] = splitSpeciesNameSuffix("CaCO3(s, calcite)"); // name: "CaCO3",   suffix: "s, calcite"
/// const auto [name, suffix] = splitSpeciesNameSuffix("MgCO3");             // name: "MgCO3",   suffix: ""
/// const auto [name, suffix] = splitSpeciesNameSuffix("Methane(aq)");       // name: "Methane", suffix: "aq"
/// const auto [name, suffix] = splitSpeciesNameSuffix("Ca++");              // name: "Ca++",    suffix: ""
/// const auto [name, suffix] = splitSpeciesNameSuffix("Fe[3+](aq)");        // name: "Fe[3+]",  suffix: "aq"
/// ~~~
auto splitSpeciesNameSuffix(std::string name) -> std::pair<std::string, std::string>;

} // namespace Reaktoro
