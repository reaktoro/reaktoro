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

namespace Reaktoro {

/// A type that describes the result of an optimisation calculation
struct OptimumResult
{
    /// The flag that indicates if the optimisation calculation converged
    bool succeeded = false;

    /// The number of iterations in the optimisation calculation
    unsigned iterations = 0;

    /// The number of evaluations of the objective function in the optimisation calculation
    unsigned num_objective_evals = 0;

    /// The convergence rate of the optimisation calculation near the solution
    double convergence_rate = 0;

    /// The final residual error of the optimisation calculation
    double error = 0;

    /// The wall time spent for the optimisation calculation (in units of s)
    double time = 0;

    /// The wall time spent for all objective evaluations (in units of s)
    double time_objective_evals = 0;

    /// The wall time spent for all contraint evaluations (in units of s)
    double time_constraint_evals = 0;

    /// The wall time spent for all linear system solutions (in units of s)
    double time_linear_systems = 0;

    /// Update this OptimumResult instance with another by addition
    auto operator+=(const OptimumResult& other) -> OptimumResult&;
};

} // namespace Reaktoro
