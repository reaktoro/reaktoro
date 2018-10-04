// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// Reaktoro includes
#include <unsupported/cpp-interpreter/Yaml.hpp>

namespace Reaktoro {

// Forward declarations
class InterpreterState;

/// Process a keyword node of type Database
auto processDatabaseNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type AqueousPhase
auto processAqueousPhaseNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type GaseousPhase
auto processGaseousPhaseNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type MineralPhase
auto processMineralPhaseNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type Minerals
auto processMineralsNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type ChemicalModel
auto processChemicalModelNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type MineralReaction
auto processMineralReactionNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type Equilibrium
auto processEquilibriumNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type EquilibriumPath
auto processEquilibriumPathNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type KineticPath
auto processKineticPathNode(InterpreterState& istate, const Node& node) -> void;

/// Process a keyword node of type Phreeqc
auto processPhreeqcNode(InterpreterState& istate, const Node& node) -> void;

} // namespace Reaktoro
