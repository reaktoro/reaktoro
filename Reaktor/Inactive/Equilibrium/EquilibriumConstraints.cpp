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

#include "EquilibriumConstraints.hpp"

// Reaktor includes
#include <Reaktor/Core/ChemicalState.hpp>

namespace Reaktor {

EquilibriumConstraints::EquilibriumConstraints(const std::function<PartialVector(const ChemicalState&)>& function, unsigned num_constraints)
: function$(function), num_constraints$(num_constraints)
{}

auto EquilibriumConstraints::numConstraints() const -> unsigned
{
    return num_constraints$;
}

auto EquilibriumConstraints::function() const -> const std::function<PartialVector(const ChemicalState&)>&
{
    return function$;
}

auto EquilibriumConstraints::operator()(const ChemicalState& state) const -> PartialVector
{
    return function$(state);
}

} /* namespace Reaktor */
