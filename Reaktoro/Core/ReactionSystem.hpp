// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// C++ includes
#include <memory>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ScalarTypes.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class ChemicalProperties;
class Reaction;

/// A class that represents a system of chemical reactions.
/// The ReactionSystem class is a collection of Reaction instances. It provides
/// convenient methods that calculates the equilibrium constants, reaction quotients,
/// and rates of the reactions.
/// @see Reaction, ChemicalSystem, ChemicalSystem
/// @ingroup Core
class ReactionSystem
{
public:
    /// Construct a default ReactionSystem instances
    ReactionSystem();

    /// Construct a ReactionSystem instance with given reactions
    ReactionSystem(const ChemicalSystem& system, const std::vector<Reaction>& reactions);

    /// Destroy this ReactionSystem instance
    virtual ~ReactionSystem();

    /// Return the number of reactions in the reaction system.
    auto numReactions() const -> unsigned;

    /// Return the index of the reaction with given name.
    auto indexReaction(std::string name) const -> Index;

    /// Return the index of the reaction with given name.
    /// It throws an exception if the element does not exist.
    auto indexReactionWithError(std::string name) const -> Index;

    /// Return the reactions in the reaction system.
    auto reactions() const -> const std::vector<Reaction>&;

    /// Return the reaction in the reaction system with given index.
    /// @param index The index of the reaction
    auto reaction(Index index) const -> const Reaction&;

    /// Return the reaction in the reaction system with given name.
    /// @param name The name of the reaction
    auto reaction(std::string name) const -> const Reaction&;

    /// Return the stoichiometric matrix of the reaction system.
    auto stoichiometricMatrix() const -> MatrixConstRef;

    /// Return the chemical system instance
    auto system() const -> const ChemicalSystem&;

    /// Calculate the equilibrium constants of the reactions.
    /// @param properties The chemical properties of the system
    auto lnEquilibriumConstants(const ChemicalProperties& properties) const -> ThermoVector;

    /// Calculate the reaction quotients of the reactions.
    /// @param properties The chemical properties of the system
    auto lnReactionQuotients(const ChemicalProperties& properties) const -> ChemicalVector;

    /// Calculate the kinetic rates of the reactions.
    /// @param properties The thermodynamic properties of the system
    auto rates(const ChemicalProperties& properties) const -> ChemicalVector;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

auto operator<<(std::ostream& out, const ReactionSystem& reactions) -> std::ostream&;


} // namespace Reaktoro
