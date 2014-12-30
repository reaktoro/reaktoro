// Reaktor is a C++ library for computational reaction modelling.
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
#include <memory>

namespace Reaktor {

// Forward declarations
class  EquilibriumProblem;
struct EquilibriumOptions;
struct EquilibriumResult;

class EquilibriumSolver
{
public:
    /// Construct a default EquilibriumSolver instance
    EquilibriumSolver();

    /// Construct a copy of an EquilibriumSolver instance
    EquilibriumSolver(const EquilibriumSolver& other);

    /// Destroy this EquilibriumSolver instance
    virtual ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver instance
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /// Find a initial guess for an equilibrium problem
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium approximation
    /// @param options The options for the equilibrium calculation
    auto approximate(const EquilibriumProblem& problem, EquilibriumResult& result) -> void;

    /// Find a initial guess for an equilibrium problem with given options
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium approximation
    /// @param options The options for the equilibrium calculation
    auto approximate(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;

    /// Solve an equilibrium problem
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium calculation
    /// @param options The options for the equilibrium calculation
    auto solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void;

    /// Solve an equilibrium problem with given options
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium calculation
    /// @param options The options for the equilibrium calculation
    auto solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
