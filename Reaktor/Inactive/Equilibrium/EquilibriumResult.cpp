/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "EquilibriumResult.hpp"

namespace Reaktor {

EquilibriumResult::EquilibriumResult()
: Optima::IPFilterResult()
{}

EquilibriumResult::EquilibriumResult(const Optima::IPFilterResult& other)
: Optima::IPFilterResult(other)
{}

EquilibriumResult operator+(const EquilibriumResult& result1, const EquilibriumResult& result2)
{
    EquilibriumResult result;
    result.converged            = result1.converged and result2.converged;
    result.num_constraint_evals = result1.num_constraint_evals + result2.num_constraint_evals;
    result.num_iterations       = result1.num_iterations + result2.num_iterations;
    result.num_objective_evals  = result1.num_objective_evals + result2.num_objective_evals;
    result.num_restorations     = result1.num_restorations + result2.num_restorations;
    return result;
}

} /* namespace Reaktor */

