// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Optimization/Hessian.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>

namespace Reaktoro {

/// A type that describes the state of an optimum solution
struct OptimumState
{
    /// The primal solution of the optimisation problem
    Vector x;

    /// The dual solution of the optimisation problem with respect to the equality constraints
    Vector y;

    /// The dual solution of the optimisation problem with respect to the lower bound constraints
    Vector z;

    /// The dual solution of the optimisation problem with respect to the upper bound constraints
    Vector w;

    /// The evaluation of the objective function at `x`
    ObjectiveResult f;
};

} // namespace Reaktoro
