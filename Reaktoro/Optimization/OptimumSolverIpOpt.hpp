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
struct OptimumResult;
struct OptimumSensitivity;
struct OptimumState;

/// The class that implements the IpOpt algorithm using an interior-point method.
/// This method is implemented based on some strategies of the optimization code IPOPT.
class OptimumSolverIpOpt : public OptimumSolverBase
{
public:
    /// Construct a default OptimumSolverIpOpt instance.
    OptimumSolverIpOpt();

    /// Construct a copy of an OptimumSolverIpOpt instance.
    OptimumSolverIpOpt(const OptimumSolverIpOpt& other);

    /// Destroy this OptimumSolverIpOpt instance.
    virtual ~OptimumSolverIpOpt();

    /// Assign an OptimumSolverIpOpt instance to this.
    auto operator=(OptimumSolverIpOpt other) -> OptimumSolverIpOpt&;

    /// Solve an optimisation problem.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Solve an optimisation problem with given options.
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation calculation
    /// @param options The options for the optimisation calculation
    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    /// Calculate the sensitivity of the optimal state with respect to a parameter *p*.
    /// @param dgdp The derivative of the gradient vector *g* with respect to the parameter *p*
    /// @param dbdp The derivative of the equality constraint vector *b* with respect to the parameter *p*
    virtual auto sensitivity(const Vector& dgdp, const Vector& dbdp) -> OptimumSensitivity;

    /// Return a clone of this instance.
    virtual auto clone() const -> OptimumSolverBase*;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
