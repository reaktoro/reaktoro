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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Optimization/Hessian.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>

namespace Reaktoro {

/// A type that describes the state of an optimum solution
struct OptimumState
{
    /// The primal solution of the optimisation problem
    VectorXd x;

    /// The dual solution of the optimisation problem with respect to the equality constraints
    VectorXd y;

    /// The dual solution of the optimisation problem with respect to the lower bound constraints
    VectorXd z;

    /// The dual solution of the optimisation problem with respect to the upper bound constraints
    VectorXd w;

    /// The evaluation of the objective function at `x`
    ObjectiveResult f;
};

} // namespace Reaktoro
