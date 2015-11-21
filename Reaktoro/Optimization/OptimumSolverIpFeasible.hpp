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

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumState;
struct OptimumResult;

/// The class that implements the IpFeasible algorithm using an interior-point method.
/// This method can be used to find an initial feasible solution.
class OptimumSolverIpFeasible
{
public:
    /// Construct a default OptimumSolverIpFeasible instance
    OptimumSolverIpFeasible();

    /// Construct a copy of an OptimumSolverIpFeasible instance
    OptimumSolverIpFeasible(const OptimumSolverIpFeasible& other);

    /// Destroy this OptimumSolverIpFeasible instance
    virtual ~OptimumSolverIpFeasible();

    /// Assign a copy of an OptimumSolverIpFeasible instance
    auto operator=(OptimumSolverIpFeasible other) -> OptimumSolverIpFeasible&;

    /// Find an initial guess for an optimisation problem
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    /// @param options The options for the optimisation calculation
    auto approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Find an initial guess for an optimisation problem with given options
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the optimisation approximation
    /// @param options The options for the optimisation calculation
    auto approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
