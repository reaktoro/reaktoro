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
#include <Reaktoro/Optimization/OptimumSolverBase.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumResult;
struct OptimumState;

/// The class that implements an optimization algorithm based on Karpov's method.
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

    /// Return the sensitivity `dx/dp` of the solution `x` with respect to a vector of parameters `p`.
    /// @param dgdp The derivatives `dg/dp` of the objective gradient `grad(f)` with respect to the parameters `p`
    /// @param dbdp The derivatives `db/dp` of the vector `b` with respect to the parameters `p`
    virtual auto dxdp(VectorConstRef dgdp, VectorConstRef dbdp) -> Vector;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase*;

  private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
