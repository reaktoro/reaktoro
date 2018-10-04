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

namespace Reaktoro {

// Forward declarations
class EquilibriumInverseProblem;
class EquilibriumProblem;
class ChemicalState;
class Partition;
struct EquilibriumOptions;
struct EquilibriumResult;

/// Equilibrate a chemical state instance.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state) -> EquilibriumResult;

/// Equilibrate a chemical state instance.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const Partition& partition) -> EquilibriumResult;

/// Equilibrate a chemical state instance with equilibrium options.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with equilibrium options.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const Partition& partition, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const EquilibriumProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(const EquilibriumProblem& problem) -> ChemicalState;

/// Equilibrate a chemical state instance with an equilibrium problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(const EquilibriumProblem& problem, const EquilibriumOptions& options) -> ChemicalState;

/// Equilibrate a chemical state instance with an equilibrium inverse problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium inverse problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem, const EquilibriumOptions& options) -> EquilibriumResult;

/// Equilibrate a chemical state instance with an equilibrium inverse problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(const EquilibriumInverseProblem& problem) -> ChemicalState;

/// Equilibrate a chemical state instance with an equilibrium inverse problem.
/// @warning This function is intended for convenience only!
/// For performance critical applications, use class EquilibriumSolver.
auto equilibrate(const EquilibriumInverseProblem& problem, const EquilibriumOptions& options) -> ChemicalState;

} // namespace Reaktoro
