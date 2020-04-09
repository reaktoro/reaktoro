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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumResult;
struct OptimumState;

/// The base class for all optimization algorithms.
class OptimumSolverBase
{
public:
    /// Pure virtual destructor
    virtual ~OptimumSolverBase() = 0;

    /// Solve an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult = 0;

    /// Solve an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    /// @param options The options for the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult = 0;

    /// Return the sensitivity `dx/dp` of the solution `x` with respect to a vector of parameters `p`.
    /// @param dgdp The derivatives `dg/dp` of the objective gradient `grad(f)` with respect to the parameters `p`
    /// @param dbdp The derivatives `db/dp` of the vector `b` with respect to the parameters `p`
    virtual auto dxdp(VectorXdConstRef dgdp, VectorXdConstRef dbdp) -> VectorXd = 0;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase* = 0;
};

} // namespace Reaktoro
