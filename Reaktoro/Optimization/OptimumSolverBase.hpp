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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumResult;
struct OptimumSensitivity;
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

    /// Calculate the sensitivity of the optimal state with respect to a parameter *p*.
    /// @param dgdp The derivative of the gradient vector *g* with respect to the parameter *p*
    /// @param dbdp The derivative of the equality constraint vector *b* with respect to the parameter *p*
    virtual auto sensitivity(const Vector& dgdp, const Vector& dbdp) -> OptimumSensitivity = 0;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase* = 0;
};

} // namespace Reaktoro
