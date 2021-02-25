// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Reaction instances for less than
auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool;

/// Compare two Reaction instances for equality
auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool;

} // namespace Reaktoro
