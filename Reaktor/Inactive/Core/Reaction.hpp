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

// Reaktor includes
#include <Reaktor/Core/Types.hpp>
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/PartialVector.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalSystem;
class ReactionEquation;

/**
 * Provides a computational representation of a chemical reaction
 *
 * The Reaction class provides a representation of a chemical reaction and operations
 * such as the calculation of equilibrium constants at given temperature and pressure points,
 * reaction quotients, and reaction rates.
 *
 * @see ReactionEquation, ReactionRate, EquilibriumConstant, ReactionSystem, ChemicalSystem
 * @ingroup Core
 */
class Reaction
{
public:
	/**
	 * Constructs a default Reaction instance
	 */
	Reaction();

	/**
     * Constructs a Reaction instance
     * @param equation The reaction equation instance
     * @param system The chemical system instance
     * @see ReactionEquation, ChemicalSystem
     */
    Reaction(const ReactionEquation& equation, const ChemicalSystem& system);

    /**
     * Constructs a copy of a Reaction instance
     */
    Reaction(const Reaction& other);

    /**
     * Destroys this instance
     */
    virtual ~Reaction();

    /**
     * Assigns a Reaction instance to this instance
     */
    auto operator=(Reaction other) -> Reaction&;

	/**
	 * Sets the equilibrium constant function of the reaction
	 *
	 * **Note:** If no equilibrium contant function is provided, a default
	 * one will be used, which is defined in terms of the chemical potentials
	 * of the participating species in the reaction.
	 *
	 * @param K The equilibrium contant function
	 */
	auto setEquilibriumConstant(const EquilibriumConstantFn& equilibrium_constant) -> void;

	/**
	 * Sets the reaction rate function of the reaction
	 *
	 * **Note:** If no reaction rate function is provided, a default
     * one will be used, which is defined as a zero rate function.
     *
     * @param rate The reaction rate instance
	 */
	auto setRate(const ReactionRateFn& rate) -> void;

	/**
	 * Gets the reaction equation of the reaction
	 */
	auto equation() const -> const ReactionEquation&;

	/**
	 * Gets the reaction rate function of the reaction
	 */
	auto rate() const -> const ReactionRateFn&;

	/**
	 * Gets the equilibrium constant function of the reaction
	 */
	auto equilibriumConstant() const -> const EquilibriumConstantFn&;

	/**
	 * Gets the indices of the reacting species in the reaction
	 */
	auto idxReactingSpecies() const -> const Indices&;

	/**
	 * Gets the stoichiometry of a species in the reaction
	 */
	auto stoichiometry(const std::string& species) const -> double;

	/**
	 * Calculates the equilibrium constant of the reaction
	 * @param T The temperature of the chemical system (in units of K)
     * @param P The pressure of the chemical system (in units of Pa)
	 */
	auto equilibriumConstant(double T, double P) const -> double;

	/**
	 * Calculates the reaction quotient of the reaction
	 *
	 * The reaction quotient @f$Q@f$ of a reaction is defined as:
	 * @f[
     *     Q=\prod_{i=1}^{N}a_{i}^{\nu_{i}},
     * @f]
     * where @f$N@f$ denotes the number of species in the chemical system,
     * @f$a_{i}@f$ the activity of the @f$i@f$-th species, and
     * @f$\nu_{i}@f$ the stoichiometry of the @f$i@f$-th species in the
     * reaction:
     * @f[
     *     0\rightleftharpoons\sum_{i=1}^{N}\nu_{i}\alpha_{i},
     * @f]
     * with @f$\alpha_{i}@f$ denoting the @f$i@f$-th species. The sign
     * convention for the stoichiometric coefficients is: *positive* for
     * products, *negative* for reactants.
     *
	 * @param a The activities of every species in the chemical system and their molar derivatives
	 */
	auto reactionQuotient(const PartialVector& a) const -> PartialScalar;

    /**
     * Calculates the kinetic rate of the reaction
     * @param T The temperature of the chemical system (in units of K)
     * @param P The pressure of the chemical system (in units of Pa)
     * @param n The molar abundance of the species in the chemical system (in units of mol)
     * @param a The activities of every species in the chemical system and their molar derivatives
     * @return The rate of the reaction and its molar derivatives
     */
    auto rate(double T, double P, const Vector& n, const PartialVector& a) const -> PartialScalar;

private:
    class Impl;

    std::unique_ptr<Impl> pimpl;
};

/**
 * Outputs the Reaction instance
 */
auto operator<<(std::ostream& out, const Reaction& reaction) -> std::ostream&;

} // namespace Reaktor
