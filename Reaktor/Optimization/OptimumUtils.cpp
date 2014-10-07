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

#include "OptimumUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>

namespace Reaktor {

auto numVariables(const OptimumProblem& problem) -> unsigned
{
    return problem.A.n_cols;
}

auto numConstraints(const OptimumProblem& problem) -> unsigned
{
    return problem.A.n_rows;
}

auto objective(const OptimumProblem& problem, const Vector& x) -> ObjectiveResult
{
    return problem.f(x);
}

auto constraint(const OptimumProblem& problem, const Vector& x) -> Vector
{
    return problem.A*x - problem.b;
}

} // namespace Reaktor
