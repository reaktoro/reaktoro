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

#pragma once

// C++ includes
#include <functional>

// Reaktor includes
#include <Reaktor/Common/PartialVector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalState;

class EquilibriumConstraints
{
public:
    /**
     * Constructs an EquilibriumConstraints instance
     * @param function The vector-valued function that defines the equilibrium constraints
     * @param num_constraints The number of equilibrium constraints
     */
    EquilibriumConstraints(const std::function<PartialVector(const ChemicalState&)>& function, unsigned num_constraints);

    /**
     * Gets the number of equilibrium constraints
     */
    auto numConstraints() const -> unsigned;

    /**
     * Gets the vector-valued function that defines the equilibrium constraints
     */
    auto function() const -> const std::function<PartialVector(const ChemicalState&)>&;

    /**
     * Calculates the residuals of the equilibrium constraints
     * @param state The chemical state of the system
     * @return The residuals of the equilibrium constraints and their partial molar derivatives
     */
    auto operator()(const ChemicalState& state) const -> PartialVector;

private:
    /// The vector-valued equilibrium constraint function
    std::function<PartialVector(const ChemicalState&)> function$;

    /// The number of constraints imposed by the vector-valued constraint function
    unsigned num_constraints$;
};

} // namespace Reaktor
