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

namespace Reaktoro {

// Forward declarations
class EquilibriumInverseProblem;
class EquilibriumProblem;
class EquilibriumState;
class Partition;
struct EquilibriumOptions;
struct EquilibriumResult;

/// Equilibrate a chemical state instance
auto equilibrate(EquilibriumState& state) -> EquilibriumResult;

/// Equilibrate a chemical state instance
auto equilibrate(EquilibriumState& state, const Partition& partition) -> EquilibriumResult;

/// Equilibrate a chemical state instance with equilibrium options
auto equilibrate(EquilibriumState& state, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with equilibrium options
auto equilibrate(EquilibriumState& state, const Partition& partition, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem
auto equilibrate(EquilibriumState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem
auto equilibrate(EquilibriumState& state, const EquilibriumProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem
auto equilibrate(const EquilibriumProblem& problem) -> EquilibriumState;

/// Equilibrate a chemical state instance with an equilibrium problem
auto equilibrate(const EquilibriumProblem& problem, const EquilibriumOptions& options) -> EquilibriumState;

/// Equilibrate a chemical state instance with an equilibrium inverse problem
auto equilibrate(EquilibriumState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium inverse problem
auto equilibrate(EquilibriumState& state, const EquilibriumInverseProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium inverse problem
auto equilibrate(const EquilibriumInverseProblem& problem) -> EquilibriumState;

/// Equilibrate a chemical state instance with an equilibrium inverse problem
auto equilibrate(const EquilibriumInverseProblem& problem, const EquilibriumOptions& options) -> EquilibriumState;

} // namespace Reaktoro
