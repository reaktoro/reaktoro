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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Optimization/Hessian.hpp>

namespace Reaktoro {

/// A type that describes the state of an optimum solution
struct OptimumState
{
    /// The primal solution of the optimisation problem
    Vector x;

    /// The dual solution of the optimisation problem with respect to the equality constraints
    Vector y;

    /// The dual solution of the optimisation problem with respect to the bound constraints
    Vector z;
};

} // namespace Reaktoro
