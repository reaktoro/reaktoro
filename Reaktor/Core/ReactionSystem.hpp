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

#pragma once

// C++ includes
#include <memory>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Core/Reaction.hpp>

namespace Reaktor {

// Forward declarations
class Reaction;

/// A struct to represent a more detailed model configuration of a ReactionSystem object.
/// @see Multiphase, Phase
/// @ingroup Core
struct ReactionSystemModel
{
    /// The function for the equilibrium constant of the reactions (in natural log).
    ThermoVectorFunction lnk;

    /// The function for the standard molar Gibbs free energies of the reactions (in units of J/mol).
    ThermoVectorFunction standard_gibbs_energy;

    /// The function for the standard molar Helmholtz free energies of the reactions (in units of J/mol).
    ThermoVectorFunction standard_helmholtz_energy;

    /// The function for the standard molar internal energies of the reactions (in units of J/mol).
    ThermoVectorFunction standard_internal_energy;

    /// The function for the standard molar enthalpies of the reactions (in units of J/mol).
    ThermoVectorFunction standard_enthalpy;

    /// The function for the standard molar entropies of the reactions (in units of J/K).
    ThermoVectorFunction standard_entropy;

    /// The function for the standard molar volumes of the reactions (in units of m3/mol).
    ThermoVectorFunction standard_volume;

    /// The function for the standard molar isobaric heat capacity of the reactions (in units of J/(mol*K)).
    ThermoVectorFunction standard_heat_capacity;

    /// The function for the molar volumes of the reactions (in units of m3/mol).
    ReactionRateVectorFunction rate;
};

/// A class that represents a system of chemical reactions.
/// The ReactionSystem class is a collection of Reaction instances. It provides
/// convenient methods that calculates the equilibrium constants, reaction quotients,
/// and rates of the reactions.
/// @see Reaction, Multiphase, ChemicalSystem
/// @ingroup Core
class ReactionSystem
{
public:
    /// Construct a default ReactionSystem instances
    ReactionSystem();

    /// Construct a ReactionSystem instance with given reactions
    ReactionSystem(const std::vector<Reaction>& reactions);

    /// Construct a ReactionSystem instance with given reactions and model configuration
    ReactionSystem(const std::vector<Reaction>& reactions, const ReactionSystemModel& model);

    /// Destroy this ReactionSystem instance
    virtual ~ReactionSystem();

    /// Return the number of reactions in the reaction system.
    auto numReactions() const -> unsigned;

    /// Return the reactions in the reaction system.
    auto reactions() const -> const std::vector<Reaction>&;

    /// Return the reaction in the reaction system with given index.
    /// @param index The index of the reaction
    auto reaction(Index index) const -> const Reaction&;

    /// Return the stoichiometric matrix of the reaction system.
    auto stoichiometricMatrix() const -> const Matrix&;

    /// Calculate the equilibrium constants of the reactions.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto lnEquilibriumConstants(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar Gibbs free energies of the reactions (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardGibbsEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar Helmholtz free energies of the reactions (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardHelmholtzEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar internal energies of the reactions (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardInternalEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar enthalpies of the reactions (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardEnthalpies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar entropies of the reactions (in units of J/K).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardEntropies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar volumes of the reactions (in units of m3/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardVolumes(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar isobaric heat capacity of the reactions (in units of J/(mol*K)).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardHeatCapacities(double T, double P) const -> ThermoVector;

    /// Calculate the kinetic rates of the reactions.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param n The molar amounts of the species (in units of mol)
    /// @param a The activities of the species and their partial derivatives
    auto rates(double T, double P, const Vector& n, const ChemicalVector& a) const -> ChemicalVector;

    /// Calculate the reaction quotients of the reactions.
    /// @param a The activities of the species and their partial derivatives
    auto lnReactionQuotients(const ChemicalVector& a) const -> ChemicalVector;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

auto operator<<(std::ostream& out, const ReactionSystem& reactions) -> std::ostream&;


} // namespace Reaktor
