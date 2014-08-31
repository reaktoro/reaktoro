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
#include <memory>

// Reaktor includes
#include <Reaktor/Common/PartialVector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class EquilibriumConstraints;
class Partitioning;

class BalanceConstraints
{
public:
    /**
     * Constructs a default BalanceConstraints instance
     * @param system The definition of the chemical system
     */
    explicit BalanceConstraints(const ChemicalSystem& system);

    /**
     * Constructs a default BalanceConstraints instance
     * @param system The definition of the chemical system
     * @param partitioning The partitioning of the species in the chemical system
     */
    BalanceConstraints(const ChemicalSystem& system, const Partitioning& partitioning);

    /**
     * Constructs a copy of a BalanceConstraints instance
     */
    BalanceConstraints(const BalanceConstraints& other);

    /**
     * Destroys the BalanceConstraints instance
     */
    virtual ~BalanceConstraints();

    /**
     * Assigns a BalanceConstraints instance to this
     */
    auto operator=(BalanceConstraints other) -> BalanceConstraints&;

    /**
     * Sets the molar abundance of the equilibrium elements
     * @param be The molar abundance of the equilibrium elements
     */
    auto setMassBalance(const Vector& be) -> void;

    /**
     * Gets the number of mass-balance and charge-balance constraints
     */
    auto numConstraints() const -> unsigned;

    /**
     * Gets the number of mass-balance constraints
     */
    auto numMassBalanceConstraints() const -> unsigned;

    /**
     * Gets the number of charge-balance constraints
     */
    auto numChargeBalanceConstraints() const -> unsigned;

    /**
     * Gets the indices of the phases whose charge-balance condition is not implicitly imposed by the mass-balance conditions
     */
    auto idxChargeImbalancedPhases() const -> const Indices&;

    /**
     * Gets the combined mass-balance and charge-balance matrix
     */
    auto balanceMatrix() const -> const Matrix&;

    /**
     * Calculates the residuals of the mass-balance and charge-balance constraints
     * @param state The chemical state of the system
     * @return The residuals of the balance constraints and its partial molar derivatives
     */
    auto operator()(const ChemicalState& state) -> PartialVector;

    /**
     * Creates an EquilibriumConstraints instance from this BalanceConstraints instance
     * @param be The molar abundance of the equilibrium elements
     * @return The EquilibriumConstraints instance as a result of the conversion
     */
    auto constraints(const Vector& be) const -> EquilibriumConstraints;

    /**
     * Converts this BalanceConstraints instance into an EquilibriumConstraints instance
     */
    operator EquilibriumConstraints() const;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
