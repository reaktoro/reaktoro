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

// Interpreter includes
#include "Yaml.hpp"

namespace Reaktoro {

// Forward declarations
class InterpreterState;

/// Process a keyword node of type MineralReactionNode
auto processMineralReactionNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type EquilibriumNode
auto processEquilibriumNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type KineticsNode
auto processKineticsNode(InterpreterState& istate, const Node& node) -> void;

} // namespace Reaktoro
