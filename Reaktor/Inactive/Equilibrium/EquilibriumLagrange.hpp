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
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Partitioning;

class EquilibriumLagrange
{
public:
    /**
     * Constructs an EquilibriumLagrange instance
     * @param system The ChemicalSystem instance that represents the chemical system
     */
    explicit EquilibriumLagrange(const ChemicalSystem& system);

    /**
     * Constructs an EquilibriumLagrange instance
     * @param system The ChemicalSystem instance that represents the chemical system
     * @param system The Partitioning instance that represents the partitioning of the chemical species
     */
    EquilibriumLagrange(const ChemicalSystem& system, const Partitioning& partitioning);

    /**
     * Constructs a copy of an EquilibriumLagrange instance
     */
    EquilibriumLagrange(const EquilibriumLagrange& other);

    /**
     * Destroy the EquilibriumLagrange instance
     */
    virtual ~EquilibriumLagrange();

    /**
     * Assign a EquilibriumLagrange instance to this instance
     */
    auto operator=(EquilibriumLagrange other) -> EquilibriumLagrange&;

    /**
     * Sets the Lagrange multipliers **y** and **z** of the equilibrium state of the chemical system
     *
     * @param state The chemical state of the system
     * @param y The Lagrange multipliers **y** of the equilibrium state of the system
     * @param z The Lagrange multipliers **z** of the equilibrium state of the system
     */
    auto setMultipliers(const ChemicalState& state, const Vector& y, const Vector& z) -> void;

    /**
     * Gets the definition of the chemical system
     * @return A ChemicalSystem instance representing the chemical system
     */
    auto system() const -> const ChemicalSystem&;

    /**
     * Gets the partitioning of the species in the chemical system
     * @return A Partitioning instance representing the partitioning of the species
     */
    auto partitioning() const -> const Partitioning&;

    /**
     * Gets the Lagrange multipliers **y** of the equilibrium state
     *
     * The vector of Lagrange multipliers **y**, whose dimension equals the number of
     * equilibrium elements, are the multipliers associated with the equality constraints
     * of the Gibbs energy minimisation problem.
     */
    auto lagrangeY() const -> const Vector&;

    /**
     * Gets the Lagrange multipliers **z** of the equilibrium state
     *
     * The vector of Lagrange multipliers **z**, whose dimension equals the number of
     * equilibrium species, are the multipliers associated with the inequality constraints
     * of the Gibbs energy minimisation problem.
     */
    auto lagrangeZ() const -> const Vector&;

    /**
     * Gets the stability indices of the equilibrium phases
     */
    auto stabilityIndicesEquilibriumPhases() const -> const Vector&;

    /**
     * Gets the stability indices of all phases
     *
     * The phases that do not contain equilibrium species, such as phases
     * that only contain kinetic species, have stability index equal to zero.
     */
    auto stabilityIndicesPhases() const -> const Vector&;

    /**
     * Gets the indices of the equilibrium phases that are stable
     *
     * These indices are relative to the set of equilibrium phases,
     * and not the set of all phases.
     *
     * @see idxStablePhases
     */
    auto idxEquilibriumStablePhases() const -> Indices;

    /**
     * Gets the indices of the equilibrium phases that are unstable
     *
     * These indices are relative to the set of equilibrium phases,
     * and not the set of all phases.
     *
     * @see idxUnstablePhases
     */
    auto idxEquilibriumUnstablePhases() const -> Indices;

    /**
     * Gets the indices of the equilibrium species that are stable
     *
     * These indices are relative to the set of equilibrium species,
     * and not the set of all species.
     *
     * @see idxStableSpecies
     */
    auto idxEquilibriumStableSpecies() const -> Indices;

    /**
     * Gets the indices of the equilibrium species that are unstable
     *
     * These indices are relative to the set of equilibrium species,
     * and not the set of all species.
     *
     * @see idxUnstableSpecies
     */
    auto idxEquilibriumUnstableSpecies() const -> Indices;

    /**
     * Gets the indices of the stable equilibrium phases in the system
     */
    auto idxStablePhases() const -> Indices;

    /**
     * Gets the indices of the unstable equilibrium phases in the system
     */
    auto idxUnstablePhases() const -> Indices;

    /**
     * Gets the indices of the stable equilibrium species in the system
     */
    auto idxStableSpecies() const -> Indices;

    /**
     * Gets the indices of the unstable equilibrium species in the system
     */
    auto idxUnstableSpecies() const -> Indices;

    /**
     * Gets the names of the stable equilibrium phases in the system
     */
    auto stablePhases() const -> std::vector<std::string>;

    /**
     * Gets the names of the unstable equilibrium phases in the system
     */
    auto unstablePhases() const -> std::vector<std::string>;

    /**
     * Gets the names of the stable equilibrium species in the system
     */
    auto stableSpecies() const -> std::vector<std::string>;

    /**
     * Gets the names of the unstable equilibrium species in the system
     */
    auto unstableSpecies() const -> std::vector<std::string>;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

} /* namespace Reaktor */
