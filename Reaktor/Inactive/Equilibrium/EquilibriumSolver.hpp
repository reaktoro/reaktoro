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

// Optima includes
#include <Optima/Optima.hpp>

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Reaktor forward declarations
class  ChemicalState;
class  ChemicalSystem;
class  EquilibriumConstraints;
class  EquilibriumLagrange;
class  Partitioning;
struct EquilibriumOptions;
struct EquilibriumParams;
struct EquilibriumResult;

class EquilibriumSolver
{
public:
    /**
     * Constructs a EquilibriumSolver instance
     * @param system The definition of the chemical system
     */
    explicit EquilibriumSolver(const ChemicalSystem& system);

    /**
     * Constructs a EquilibriumSolver instance
     * @param system The definition of the chemical system
     * @param partitioning The partitioning of the species in the system
     */
    EquilibriumSolver(const ChemicalSystem& system, const Partitioning& partitioning);

    /**
     * Constructs a copy of a EquilibriumSolver instance
     */
    EquilibriumSolver(const EquilibriumSolver& other);

    /**
     * Destroys the EquilibriumSolver instance
     */
    virtual ~EquilibriumSolver();

    /**
     * Assigns a EquilibriumSolver instance to this instance
     */
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /**
     * Sets the options of the chemical equilibrium solver
     */
    auto setOptions(const EquilibriumOptions& options) -> void;

    /**
     * Sets the parameters of the chemical equilibrium solver
     */
    auto setParams(const EquilibriumParams& params) -> void;

    /**
     * Sets the scaling for the chemical equilibrium calculation
     */
    auto setScaling(const Optima::Scaling& scaling) -> void;

    /**
     * Sets the scaling for the chemical equilibrium calculation
     *
     * The scaling of the variables is set to the molar abundance of the equilibrium
     * species in the chemical state.
     *
     * @param state The chemical state of the system
     */
    auto setScaling(const ChemicalState& state) -> void;

    /**
     * Minimises the Gibbs free energy of the system subject to general equilibrium constraints
     *
     * @param state The equilibrium state of chemical system
     * @param constraints The equilibrium constraints of the problem
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult;

    /**
     * Minimises the Gibbs free energy of the system subject to mass-balance constraints
     *
     * @param state The equilibrium state of chemical system
     * @param be The molar abundance of the equilibrium elements
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto minimise(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult;

    /**
     * Refines the equilibrium state of the chemical system by removing the unstable phases
     *
     * The refinement occurs by performing an equilibrium calculation using a law of mass action (LMA)
     * approach with the mass and charge balance equations as the equilibrium constraints. For this,
     * the molar abundance of the elements in the equilibrium species is retrieved from the chemical
     * state of the system.
     *
     * @param state The equilibrium state of chemical system
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto refine(ChemicalState& state, EquilibriumLagrange& lagrange) -> EquilibriumResult;

    /**
     * Refines the equilibrium state of the chemical system by removing the unstable phases
     *
     * The refinement occurs by performing an equilibrium calculation using a Gibbs energy minimisation (GEM)
     * approach with given equilibrium constraints.
     *
     * @param state The equilibrium state of chemical system
     * @param constraints The equilibrium constraints of the problem
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto refine(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult;

    /**
     * Solves an equilibrium problem by minimising the Gibbs energy of the system and performing an equilibrium refinement
     *
     * @param state The equilibrium state of chemical system
     * @param problem The equilibrium problem where the the equilibrium constraints have been specified
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto solve(ChemicalState& state, EquilibriumLagrange& lagrange, const EquilibriumConstraints& constraints) -> EquilibriumResult;

    /**
     * Solves an equilibrium problem by minimising the Gibbs energy of the system and performing an equilibrium refinement
     *
     * @param state The equilibrium state of chemical system
     * @param be The molar abundance of the equilibrium elements
     * @return The result of the equilibrium calculation as a EquilibriumResult instance
     */
    auto solve(ChemicalState& state, EquilibriumLagrange& lagrange, const Vector& be) -> EquilibriumResult;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

} /* namespace Reaktor */
