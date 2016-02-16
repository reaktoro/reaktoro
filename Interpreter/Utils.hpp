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
auto collectCompoundsInEquilibriumNode(std::set<std::string>& set, const Node& node) -> void;

/// Collect all compound names in a Speciation node.
auto collectCompoundsInSpeciationNode(std::set<std::string>& set, const Node& node) -> void;

/// Collect all compound names in all Equilibrium nodes.
auto collectCompoundsInAllEquilibriumNodes(std::set<std::string>& set, const Node& root) -> void;

/// Collect all compound names in all Speciation nodes.
auto collectCompoundsInAllSpeciationNodes(std::set<std::string>& set, const Node& root) -> void;

} // namespace Reaktoro
