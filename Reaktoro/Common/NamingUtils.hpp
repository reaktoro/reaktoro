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

#pragma once

// C++ includes
#include <string>
#include <vector>

namespace Reaktoro {

/// Return a collection of alternative names to the given species name.
/// @param name The name of the species.
auto alternativeWaterNames() -> std::vector<std::string>&;

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

/// Return the name of a charged species without charge suffix.
/// For example, `Ca++` results in `Ca`, and `HCO3-` in `HCO3`.
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

} // namespace Reaktoro
