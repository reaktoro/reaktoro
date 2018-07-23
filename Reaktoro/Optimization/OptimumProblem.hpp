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
#include <functional>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Optimization/Hessian.hpp>

namespace Reaktoro {

/// A type that describes the result of the evaluation of an objective function
struct ObjectiveResult
{
    /// The value of the objective function evaluated at `x`.
    double val = 0.0;

    /// The gradient of the objective function evaluated at `x`.
    Vector grad;

    /// The Hessian of the objective function evaluated at `x`.
    Hessian hessian;
};

/// A type that describes the functional signature of an objective function.
/// @param x The vector of primal variables
/// @return The objective function evaluated at `x`
using ObjectiveFunction = std::function<ObjectiveResult(const Vector& x)>;

/// A type that describes the non-linear constrained optimisation problem
struct OptimumProblem
{
    /// The objective function.
    ObjectiveFunction objective;

    /// The number of primal variables `x`
    Index n;

    /// The coefficient vector of a linear programming problem `min tr(c)*x subject to A*x = b`.
    Vector c;

    /// The coefficient matrix of the linear equality constraint `A*x = b`.
    Matrix A;

    /// The right-hand side vector of the linear equality constraint `A*x = b`.
    Vector b;

    /// The coefficient matrix of the linear inequality constraint `Ai*x = bi`.
    Matrix Ai;

    /// The right-hand side vector of the linear equality constraint `Ai*x = bi`.
    Vector bi;

    /// The lower bound of the primal variables `x`.
    Vector l;

    /// The upper bound of the primal variables `x`.
    Vector u;
};

/// Returns true if the evaluation of a objective function has finite value and gradient.
auto isfinite(const ObjectiveResult& f) -> bool;

} // namespace Reaktoro
