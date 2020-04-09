// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
struct OptimumSensitivity;
struct OptimumState;

/// The class that implements the ActNewton algorithm using an active-set strategy.
class OptimumSolverActNewton : public OptimumSolverBase
{
public:
    /// Construct a default OptimumSolverActNewton instance.
    OptimumSolverActNewton();

    /// Construct a copy of an OptimumSolverActNewton instance.
    OptimumSolverActNewton(const OptimumSolverActNewton& other);

    /// Destroy this OptimumSolverActNewton instance.
    virtual ~OptimumSolverActNewton();

    /// Assign an OptimumSolverActNewton instance to this.
    auto operator=(OptimumSolverActNewton other) -> OptimumSolverActNewton&;

    /// Solve an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Solve an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    /// @param options The options for the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    /// Return the sensitivity `dx/dp` of the solution `x` with respect to a vector of parameters `p`.
    /// @param dgdp The derivatives `dg/dp` of the objective gradient `grad(f)` with respect to the parameters `p`
    /// @param dbdp The derivatives `db/dp` of the vector `b` with respect to the parameters `p`
    virtual auto dxdp(VectorXdConstRef dgdp, VectorXdConstRef dbdp) -> VectorXd;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase*;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
