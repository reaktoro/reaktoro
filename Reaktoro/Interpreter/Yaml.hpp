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

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

namespace Reaktoro {

/// The alias type to YAML::Node.
using Node = YAML::Node;

/// The type used to represent a node processor function.
using ProcessFunction = std::function<void(const Node&)>;

/// Return a string representation of a YAML node.
auto str(const Node& node) -> std::string;

/// Return the left-hand side node of a YAML map node.
auto keynode(const Node& node) -> Node;

/// Return the right-hand side node of a YAML map node.
auto valnode(const Node& node) -> Node;

/// Return the first word before space on the left-hand side of a YAML map node.
auto keyword(const Node& node) -> std::string;

/// Return the second word after space on the left-hand side of a YAML map node.
auto identifier(const Node& node) -> std::string;

/// Return a joined string with a node string representation.
auto operator+(std::string str, const Node& node) -> std::string;

/// Return a joined string with a node string representation.
auto operator+(const Node& node, std::string str) -> std::string;

} // namespace Reaktoro
