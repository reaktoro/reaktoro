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

// Reaktoro includes
#include <unsupported/cpp-interpreter/Keywords.hpp>
#include <unsupported/cpp-interpreter/Yaml.hpp>

namespace Reaktoro {
namespace kwd {

/// Convert a YAML node to an instance of std::string.
auto operator>>(const Node& node, std::string& x) -> void;

/// Convert a YAML node to an instance of ValueUnits.
auto operator>>(const Node& node, ValueUnits& x) -> void;

/// Convert a YAML node to an instance of EntityValueUnits.
auto operator>>(const Node& node, EntityValueUnits& x) -> void;

/// Convert a YAML node to an instance of ValueUnitsEntity.
auto operator>>(const Node& node, ValueUnitsEntity& x) -> void;

/// Convert a YAML node to an instance of Recipe.
auto operator>>(const Node& node, Recipe& x) -> void;

/// Convert a YAML node to an instance of pH.
auto operator>>(const Node& node, pH& x) -> void;

/// Convert a YAML node to an instance of SpeciesAmount.
auto operator>>(const Node& node, SpeciesAmount& x) -> void;

/// Convert a YAML node to an instance of SpeciesActivity.
auto operator>>(const Node& node, SpeciesActivity& x) -> void;

/// Convert a YAML node to an instance of SpeciesFugacity.
auto operator>>(const Node& node, SpeciesFugacity& x) -> void;

/// Convert a YAML node to an instance of PhaseAmount.
auto operator>>(const Node& node, PhaseAmount& x) -> void;

/// Convert a YAML node to an instance of PhaseVolume.
auto operator>>(const Node& node, PhaseVolume& x) -> void;

/// Convert a YAML node to an instance of PlotKeyword.
auto operator>>(const Node& node, Plot& x) -> void;

/// Convert a YAML node to an instance of EquilibriumProblem.
auto operator>>(const Node& node, EquilibriumProblem& x) -> void;

/// Convert a YAML node to an instance of EquilibriumPath.
auto operator>>(const Node& node, EquilibriumPath& x) -> void;

/// Convert a YAML node to an instance of KineticPath.
auto operator>>(const Node& node, KineticPath& x) -> void;

/// Convert a YAML node to an instance of MineralReaction.
auto operator>>(const Node& node, MineralReaction& x) -> void;

/// Convert a YAML node to an instance of SpeciationInput.
auto operator>>(const Node& node, Concentrations& x) -> void;

/// Convert a YAML node to an instance of SpeciationProblem.
auto operator>>(const Node& node, SpeciationProblem& x) -> void;

/// Convert a YAML node to an instance of PhreeqcKeyword.
auto operator>>(const Node& node, PhreeqcKeyword& x) -> void;

} // namespace kwd
} // namespace Reaktoro
