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
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

/// A type that describes the result of the evaluation of a objective function
/// @see ObjectiveFunction
struct ObjectiveResult
{
    /// The result of the evaluation of the objective function
    double func = 0.0;

    /// The gradient result of the evaluation of the objective function
    Vector grad;

    /// The Hessian result of the evaluation of the objective function
    Matrix hessian;
};

/// A type that describes the result of the evaluation of a constraint function
/// @see ConstraintFunction
struct ConstraintResult
{
    /// The result of the evaluation of the constraint function
    Vector func;

    /// The gradient result of the evaluation of the constraint function
    Matrix grad;
};

/// A type that describes the functional signature of an objective function
/// @see ObjectiveResult
typedef std::function<
    ObjectiveResult(const Vector&)>
        ObjectiveFunction;

/// A type that describes the functional signature of a constraint function
/// @see ConstraintResult
typedef std::function<
    ConstraintResult(const Vector&)>
        ConstraintFunction;

/// A type that describes the definition of an optimisation problem
class OptimumProblem
{
public:
    /// Construct a OptimumProblem instance
    /// @param num_variables The number of primal variables
    /// @param num_constraints The number of equality constraints
    OptimumProblem(unsigned num_variables, unsigned num_constraints);

    /// Set the objective function of the optimisation problem
    /// @param f The objective function of the optimisation problem
    auto setObjective(const ObjectiveFunction& objective) -> void;

    /// Set the equality constraint function of the optimisation problem
    /// @param A The coefficient matrix of the equality constraints
    /// @param b The right-hand side vector of the equality constraints
    auto setConstraint(const ConstraintFunction& constraint) -> void;

    /// Set the lower bounds of the optimisation problem
    /// @param l The lower bounds of the primal variables
    auto setLowerBounds(const Vector& lower) -> void;

    /// Set the lower bounds of the optimisation problem
    /// @param l The lower bounds of the primal variables
    auto setLowerBounds(double lower) -> void;

    /// Set the lower bounds of the optimisation problem
    /// @param u The upper bounds of the primal variables
    auto setUpperBounds(const Vector& upper) -> void;

    /// Set the lower bounds of the optimisation problem
    /// @param u The upper bounds of the primal variables
    auto setUpperBounds(double upper) -> void;

    /// Get the number of variables in the optimisation problem
    auto numVariables() const -> unsigned;

    /// Get the number of equality constraints in the optimisation problem
    auto numConstraints() const -> unsigned;

    /// Get the objective function of the optimisation problem
    auto objective() const -> const ObjectiveFunction&;

    /// Get the equality constraint function of the optimisation problem
    auto constraint() const -> const ConstraintFunction&;

    /// Get the lower bounds of the optimisation problem
    auto lowerBounds() const -> const Vector&;

    /// Get the upper bounds of the optimisation problem
    auto upperBounds() const -> const Vector&;

private:
    /// The number of primal variables
    unsigned n;

    /// The number of equality constraints
    unsigned m;

    /// The objective function to be minimized in the optimisation problem
    ObjectiveFunction f;

    /// The coefficient matrix of the equality constraints in the optimisation problem
    ConstraintFunction h;

    /// The lower bound vector of the inequality constraints in the optimisation problem
    Vector l;

    /// The upper bound vector of the inequality constraints in the optimisation problem
    Vector u;
};

} // namespace Reaktor
