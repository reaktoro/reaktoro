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

#pragma once

// C++ includes
#include <istream>
#include <string>
#include <set>

// Interpreter includes
#include "Yaml.hpp"

namespace Reaktoro {

/// Return a preprocessed input script that conforms with YAML rules.
auto preprocess(std::string script) -> std::string;

/// Return a preprocessed input script that conforms with YAML rules.
auto preprocess(std::istream& stream) -> std::string;

/// Collect all compound names in an Equilibrium node.
/// @param node The YAML node describing an Equilibrium block.
auto collectCompoundsInEquilibriumNode(std::set<std::string>& set, const Node& node) -> void;

/// Collect all compound names in a Speciation node.
/// @param node The YAML node describing a Speciation block.
auto collectCompoundsInSpeciationNode(std::set<std::string>& set, const Node& node) -> void;

/// Collect all compound names in in a input script file.
/// @param root The YAML root node of the script file.
auto collectCompounds(const Node& root) -> std::vector<std::string>;

/// Return all elements that compose the compounds in in a input script file.
/// A Database argument is needed because some compounds can be species names,
/// e.g., Calcite, Quartz, whose constituent elements are not immediately determined
/// from their names. This would not be the case for `CaCO3` and `SiO2`.
/// @param compounds The names of the compounds found in a script file.
/// @param database The database used to identify the elements of compounds that are species names found in the database.
auto identifyElements(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

/// Return a list of compounds names that are present as gaseous species in a database.
auto filterGaseousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

/// Return a list of compounds names that are present as mineral species in a database.
auto filterMineralSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

/// Return true if the input script has `Speciation` keyword.
/// This method is used to determine if all minerals in a database that could be formed from
/// a list of compounds should be returned.
/// @see speciate
auto hasSpeciation(const Node& root) -> bool;

/// Return a list of aqueous, gaseous, and mineral species that can be formed out of a given list of compounds.
/// The returned aqueous species are all those in the database whose elements are those found in the list of compounds.
/// The returned gaseous and mineral species are all those compounds whose names can be found in the database.
/// Note, however, that if the input script contains Speciation definition, then all minerals in the database that
/// could be formed out of the list of compounds is returned.
auto speciate(const Node& root, const Database& database,
    std::vector<std::string>& aqueous_species,
    std::vector<std::string>& gaseous_species,
    std::vector<std::string>& mineral_species) -> void;

/// Return a list of aqueous species that can be formed out of a given list of compounds.
auto identifyAqueousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

/// Return a list of gaseous species that can be formed out of a given list of compounds.
auto identifyGaseousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

/// Return a list of mineral species that can be formed out of a given list of compounds.
auto identifyMineralSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>;

} // namespace Reaktoro
