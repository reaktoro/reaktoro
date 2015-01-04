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
struct OptimumOptions;
class  OptimumProblem;
struct OptimumState;
struct OptimumResult;

class OptimumSolverIpfeasible
{
public:
    /// Construct a default OptimumSolverIpfeasible instance
    OptimumSolverIpfeasible();

    /// Construct a copy of an OptimumSolverIpfeasible instance
    OptimumSolverIpfeasible(const OptimumSolverIpfeasible& other);

    /// Destroy this OptimumSolverIpfeasible instance
    virtual ~OptimumSolverIpfeasible();

    /// Assign a copy of an OptimumSolverIpfeasible instance
    auto operator=(OptimumSolverIpfeasible other) -> OptimumSolverIpfeasible&;

    /// Find a initial guess for an optimisation problem
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param options The options for the optimisation calculation
    auto approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    /// Find a initial guess for an optimisation problem with given options
    /// @param problem The definition of the optimisation problem
    /// @param state[in,out] The initial guess and the final state of the equilibrium approximation
    /// @param options The options for the optimisation calculation
    auto approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
