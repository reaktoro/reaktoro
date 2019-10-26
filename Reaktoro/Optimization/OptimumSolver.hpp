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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumResult;
struct OptimumSensitivity;
struct OptimumState;
enum class OptimumMethod;

/// The friendly interface to all optimisation algorithms.
class OptimumSolver
{
public:
    /// Construct a default OptimumSolver instance.
    OptimumSolver();

    /// Construct an OptimumSolver instance with given method.
    OptimumSolver(OptimumMethod method);

    /// Construct a copy of an OptimumSolver instance.
    OptimumSolver(const OptimumSolver& other);

    /// Destroy this OptimumSolver instance.
    virtual ~OptimumSolver();

    /// Assign a copy of an OptimumSolver instance.
    auto operator=(OptimumSolver other) -> OptimumSolver&;

    /// Set the optimisation method.
    auto setMethod(OptimumMethod method) -> void;

    /// Find an initial guess for an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    auto approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Find an initial guess for an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    /// @param options The options for the optimisation calculation
    auto approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    /// Solve an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Solve an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    /// @param options The options for the optimisation calculation
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    /// Return the sensitivities `dx/dp`, `dy/dp`, `dz/dp` of the solution `(x,y,z)` with respect to a vector of parameters `p`.
    /// @param dgdp The derivatives `dg/dp` of the objective gradient `grad(f)` with respect to the parameters `p`
    /// @param dbdp The derivatives `db/dp` of the vector `b` with respect to the parameters `p`
    /// @param[out] dxdp The derivatives `dx/dp`
    /// @param[out] dydp The derivatives `dy/dp`
    /// @param[out] dzdp The derivatives `dz/dp`
    auto sensitivities(const Matrix& dgdp, const Matrix& dbdp, MatrixRef dxdp, MatrixRef dydp, MatrixRef dzdp) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
