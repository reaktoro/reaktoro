// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;

/// The function type for calculation of equilibrium constant of reaction (in natural log).
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
/// @return The ln equilibrium constant of the reaction (@eq{\ln{K}})
using EquilibriumConstantFn = std::function<real(real T, real P)>;

/// The function type for calculation of reaction rates (in mol/s).
/// @param props The chemical properties of the chemical system
/// @return The rate of the reaction (in mol/s)
/// @see Reaction
/// @ingroup Core
using ReactionRateFn = std::function<real(const ChemicalProps& props)>;

/// A class to represent a reaction and its attributes.
/// The Reaction class provides a representation of a chemical reaction and
/// operations such as the calculation of equilibrium constants at given
/// temperature and pressure points, reaction quotients, and reaction rates.
/// @see ReactionRate, EquilibriumConstant
/// @ingroup Core
class Reaction
{
public:
    /// Construct a default Reaction instance
    Reaction();

    /// Construct a Reaction instance from a ReactionEquation instance
    Reaction(const ReactionEquation& equation, const ChemicalSystem& system);

    /// Return a deep copy of this Reaction object.
    auto clone() const -> Reaction;

    /// Return a duplicate of this Reaction object with new reaction name.
    auto withName(String name) const -> Reaction;

    /// Return a duplicate of this Reaction object with new reaction equation.
    auto withEquation(const ReactionEquation& equation) const -> Reaction;

    /// Return a duplicate of this Reaction object with new equilibrium constant function.
    auto withEquilibriumConstantFn(const EquilibriumConstantFn& fn) const -> Reaction;

    /// Return a duplicate of this Reaction object with new reaction rate function.
    auto withRateFn(const ReactionRateFn& fn) const -> Reaction;

    /// Return the name of the reaction.
    auto name() const -> String;

    /// Return the equation of the reaction.
    auto equation() const -> const ReactionEquation&;

    /// Return the equilibrium constant function of the reaction.
    auto equilibriumConstantFn() const -> const EquilibriumConstantFn&;

    /// Return the rate function of the reaction.
    auto rateFn() const -> const ReactionRateFn&;

    /// Calculate the equilibrium constant of the reaction (in natural log).
    /// @param properties The chemical properties of the system
    auto lnEquilibriumConstant(const ChemicalProps& properties) const -> real;

    /// Calculate the reaction quotient of the reaction (in natural log scale).
    /// The reaction quotient of a reaction is defined as:
    /// @f[\ln Q=\sum_{i=1}^{N}\nu_{i}\ln a_{i},@f]
    /// where @f$N@f$ denotes the number of species in the multiphase system,
    /// @f$a_{i}@f$ the activity of the @f$i@f$-th species, and
    /// @f$\nu_{i}@f$ the stoichiometry of the @f$i@f$-th species in the reaction:
    /// @f[0\rightleftharpoons\sum_{i=1}^{N}\nu_{i}\alpha_{i},@f]
    /// with @f$\alpha_{i}@f$ denoting the @f$i@f$-th species. The sign
    /// convention for the stoichiometric coefficients is: *positive* for
    /// products, *negative* for reactants.
    /// @param properties The chemical properties of the system
    auto lnReactionQuotient(const ChemicalProps& properties) const -> real;

    /// Calculate the equilibrium index of the reaction as @f$\ln(Q/K)@f$.
    auto lnEquilibriumIndex(const ChemicalProps& properties) const -> real;

    /// Return the stoichiometry of a species in the reaction equation.
    /// @param species The name of the species.
    auto stoichiometry(std::string species) const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Reaction instances for less than
auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool;

/// Compare two Reaction instances for equality
auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool;

} // namespace Reaktoro
