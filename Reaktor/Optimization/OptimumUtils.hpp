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
#include <functional>

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

/// A type that describes the result of the evaluation of a objective function
/// @see ObjectiveFunction
struct ObjectiveResult
{
    /// The result of the evaluation of the objective function
    double func;

    /// The gradient result of the evaluation of the objective function
    Vector grad;

    /// The Hessian result of the evaluation of the objective function
    Matrix hessian;
};

/// A type that describes the functional signature of an objective function
/// @see ObjectiveResult
typedef std::function<
    ObjectiveResult(const Vector&)>
        ObjectiveFunction;

/// A type that describes the solution of an optimization problem
struct OptimumSolution
{
    /// The primal solution of the optimization problem
    Vector x;

    /// The dual solution of the optimization problem with respect to the equality constraints
    Vector y;

    /// The dual solution of the optimization problem with respect to the lower bound constraints
    Vector zl;

    /// The dual solution of the optimization problem with respect to the upper bound constraints
    Vector zu;
};

/// A type that describes the result of an optimization calculation
struct OptimumStatistics
{
    /// The flag that indicates if the optimization calculation converged
    bool converged;

    /// The number of iterations in the optimization calculation
    unsigned num_iterations;

    /// The number of evaluations of the objective function in the optimization calculation
    unsigned num_objective_evals;

    /// The convergence rate of the optimization calculation near the solution
    double convergence_rate;

    /// The final residual error of the optimization calculation
    double error;
};

/// A type that contains important internal state for the optimization calculation.
/// @warning The internal state of this class must not be changed outside the optimization algorithms.
struct OptimumInternal
{

};

/// A type that describes the result of an optimization calculation
struct OptimumResult
{
    /// The solution of the optimization calculation
    OptimumSolution solution;

    /// The statistics of the optimization calculation
    OptimumStatistics statistics;

    /// The internal state of the optimization calculation
    OptimumInternal internal;
};

/// A type that describes the definition of an optimization problem
struct OptimumProblem
{
    /// The objective function to be minimized in the optimization problem
    ObjectiveFunction f;

    /// The coefficient matrix of the equality constraints in the optimization problem
    Matrix A;

    /// The right-hand side vector of the equality constraints in the optimization problem
    Vector b;

    /// The lower bound vector of the inequality constraints in the optimization problem
    Vector l;

    /// The upper bound vector of the inequality constraints in the optimization problem
    Vector u;
};

/// A type that describes
struct OptimumOptions
{
    /// The residual tolerance in the optimization calculation
    double tolerance = 1.0e-8;

    /// The interior-point method perturbation in the optimization calculation
    double mu = 1.0e-8;
};

/// Get the number of variables in an optimization problem
/// @param problem The optimization problem instance
auto numVariables(const OptimumProblem& problem) -> unsigned;

/// Get the number of equality constraints in an optimization problem
/// @param problem The optimization problem instance
auto numConstraints(const OptimumProblem& problem) -> unsigned;

/// Evaluate the objective function of a optimization problem
/// @param problem The optimization problem instance
/// @param x The point where the objective function is evaluated
auto objective(const OptimumProblem& problem, const Vector& x) -> ObjectiveResult;

/// Evaluate the constraint function of a optimization problem
/// @param problem The optimization problem instance
/// @param x The point where the constraint function is evaluated
auto constraint(const OptimumProblem& problem, const Vector& x) -> Vector;

} // namespace Reaktor
