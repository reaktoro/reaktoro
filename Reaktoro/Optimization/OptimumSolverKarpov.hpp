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

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumSolverBase.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumState;
struct OptimumResult;

class OptimumSolverKarpov : public OptimumSolverBase
{
public:
    /// Construct a default OptimumSolverKarpov instance
    OptimumSolverKarpov();

    /// Construct a copy of an OptimumSolverKarpov instance
    OptimumSolverKarpov(const OptimumSolverKarpov& other);

    /// Destroy this OptimumSolverKarpov instance
    virtual ~OptimumSolverKarpov();

    /// Assign a copy of an OptimumSolverKarpov instance
    auto operator=(OptimumSolverKarpov other) -> OptimumSolverKarpov&;

    /// Solve the linear optimisation problem by finding a feasible point and then applying a simplex algorithm.
    /// @param problem The definition of the linear optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Solve the linear optimisation problem by finding a feasible point and then applying a simplex algorithm.
    /// @param problem The definition of the linear optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    /// @param options The options for the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    virtual auto clone() const -> OptimumSolverBase*;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
