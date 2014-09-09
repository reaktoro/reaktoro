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
#include <iostream>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/VectorResult.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Reaction.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalState;
class ChemicalSystem;

/**
 * Provides a computational representation of a system of chemical reactions
 *
 * The ReactionSystem class is a collection of Reaction instances. It provides
 * convenient methods that calculates the equilibrium constants, reaction quotients,
 * and rates of the reactions.
 *
 * @see Reaction, ChemicalSystem
 * @ingroup Core
 */
class ReactionSystem
{
public:
	/**
	 * Constructs a default ReactionSystem instances
	 */
	ReactionSystem();

	/**
	 * Constructs a ReactionSystem instance with given reactions
	 * @param reactions The chemical reactions instances
	 * @param system The chemical system instance
	 */
	ReactionSystem(const std::vector<Reaction>& reactions, const ChemicalSystem& system);

	/**
	 * Gets the number of reactions in the reaction system
	 */
	auto numReactions() const -> unsigned;

	/**
	 * Gets the reactions in the reaction system
	 */
	auto reactions() const -> const std::vector<Reaction>&;

	/**
	 * Gets the stoichiometric matrix of the reaction system
	 */
	auto stoichiometricMatrix() const -> const Matrix&;

	/**
     * Gets the mapping from reaction to phase
     *
     * The mapping *reaction-to-phases* can be used to determine
     * the indices of the phases involved in each reaction.
     *
     * **Usage**
     *
     * Below we show how to output the involved phases in each chemical reaction.
     * For this, we let `system` denote an instance of class ChemicalSystem, and
     * `reactions` an instance of class ReactionSystem.
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * // The mapping reaction-to-phase of the reaction system
     * const Indices& map = reactions.mapReactionToPhase();
     *
     * // Iterate over all reactions in the reaction system
     * for(Index i = 0; i < reactions.numReactions(); ++i)
     * {
     *     std::cout << "The phases involved in reaction " << reactions.reaction(i) << " are: ";
     *     for(Index j : map[i])
     *         std::cout << system.phase(j).name() << " ";
     *     std::cout << std::endl;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * @return The indices of the involved phases in each reaction
     */
	auto mapReactionToPhase() const -> const std::vector<Indices>&;

	/**
	 * Gets the mapping from reaction to species
	 *
     * The mapping *reaction-to-species* can be used to determine
     * the indices of the species involved in each reaction.
     *
     * **Usage**
     *
     * Below we show how to output the involved species in each chemical reaction.
     * For this, we let `system` denote an instance of class ChemicalSystem, and
     * `reactions` an instance of class ReactionSystem.
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * // The mapping reaction-to-species of the reaction system
     * const Indices& map = reactions.mapReactionToSpecies();
     *
     * // Iterate over all reactions in the reaction system
     * for(Index i = 0; i < reactions.numReactions(); ++i)
     * {
     *     std::cout << "The species involved in reaction " << reactions.reaction(i) << " are: ";
     *     for(Index j : map[i])
     *         std::cout << system.species(j).name() << " ";
     *     std::cout << std::endl;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * @return The indices of the involved species in each reaction
     */
	auto mapReactionToSpecies() const -> const std::vector<Indices>&;

	/**
     * Gets the mapping from species to reactions
     *
     * The mapping *species-to-reactions* can be used to determine
     * the indices of the reactions which each species is involved.
     *
     * **Usage**
     *
     * Below we show how to output the chemical reactions where each species is involved.
     * For this, we let `system` denote an instance of class ChemicalSystem, and
     * `reactions` an instance of class ReactionSystem.
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * // The mapping species-to-reactions of the reaction system
     * const Indices& map = reactions.mapSpeciesToReactions();
     *
     * // Iterate over all species in the system
     * for(Index i = 0; i < system.numSpecies(); ++i)
     * {
     *     std::cout << "The species " << system.species(i).name() << " participates in the following reactions: " << std::endl;
     *     for(Index j : map[i])
     *         std::cout << reactions.reaction(j) << std::endl;
     *     std::cout << std::endl;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * @return The indices of the involved species in each reaction
     */
	auto mapSpeciesToReaction() const -> const std::vector<Indices>&;

	/**
	 * Calculates the equilibrium constant of the reactions in reaction system
	 * @param T The temperature of the chemical system (in units of K)
	 * @param P The pressure of the chemical system (in units of Pa)
	 * @see Reaction::equilibriumConstant
	 */
	auto equilibriumConstants(double T, double P) const -> Vector;

	/**
	 * Calculates the reaction quotients of the reactions in the reaction system
	 * @param a The activities of the species in the chemical system and its molar derivatives
	 * @see Reaction::reactionQuotient
	 */
	auto reactionQuotients(const VectorResult& a) const -> VectorResult;

	/**
     * Calculates the kinetic rates of the reactions in the reaction system
     * @param T The temperature of the chemical system (in units of K)
     * @param P The pressure of the chemical system (in units of Pa)
     * @param n The molar abundance of the species in the chemical system (in units of mol)
     * @param a The activities of every species in the chemical system and their molar derivatives
     * @return The rates of the reactions and their molar derivatives
     * @see Reaction::rate
     */
	auto rates(double T, double P, const Vector& n, const VectorResult& a) const -> VectorResult;

	/**
     * Calculates the kinetic rates of the reactions in the reaction system
     * @param state The state of the chemical system
     * @param a The activities of every species in the chemical system and their molar derivatives
     * @return The rates of the reactions and their molar derivatives
     * @see Reaction::rate, ChemicalState
     */
	auto rates(const ChemicalState& state, const VectorResult& a) const -> VectorResult;

private:
	/// The reactions that compose the reaction system
	std::vector<Reaction> m_reactions;

	/// The stoichiometric matrix of the system of reactions w.r.t. to all species in the chemical system (not only the reactive ones)
	Matrix stoichiometric_matrix;

	/// The mapping from a chemical reaction to the indices of the phases that participates in it
	std::vector<Indices> reaction_phase_map;

	/// The mapping from a chemical reaction to the indices of the chemical species that participates in it
	std::vector<Indices> reaction_species_map;

	/// The mapping from a chemical species to the indices of the chemical reactions where it participates in
	std::vector<Indices> species_reaction_map;

private:
	auto initialiseStoichiometricMatrix(const ChemicalSystem& system) -> void;

	auto initialiseReactionPhaseMap(const ChemicalSystem& system) -> void;

	auto initialiseReactionSpeciesMap() -> void;

	auto initialiseSpeciesReactionMap(const ChemicalSystem& system) -> void;
};

auto operator<<(std::ostream& out, const ReactionSystem& reactions) -> std::ostream&;

} // namespace Reaktor
