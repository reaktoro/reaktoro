// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ReactionProps.hpp>
#include <Reaktoro/Core/ReactionRateModel.hpp>

namespace Reaktoro {

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

    /// Return a duplicate of this Reaction object with new reaction rate model.
    auto withRateModel(const ReactionRateModel& model) const -> Reaction;

    /// Return the name of the reaction.
    auto name() const -> String;

    /// Return the equation of the reaction.
    auto equation() const -> const ReactionEquation&;

    /// Return the rate model of the reaction.
    auto rateModel() const -> const ReactionRateModel&;

    /// Calculate the complete set of thermodynamic properties of the reaction.
    /// @param T The temperature for the calculation (in K)
    /// @param P The pressure for the calculation (in Pa)
    auto props(real T, real P) const -> ReactionProps;

    /// Calculate the complete set of thermodynamic properties of the reaction.
    /// @param T The temperature for the calculation
    /// @param unitT The temperature unit for the calculation
    /// @param P The pressure for the calculation
    /// @param unitP The pressure unit for the calculation
    auto props(real T, Chars unitT, real P, Chars unitP) const -> ReactionProps;

    /// Calculate the rate of the reaction for given chemical state.
    /// @attention Make sure `state` has updated chemical properties before
    /// using it to evaluate reaction rates. If this is not the case,
    /// ensure `state.props().update(state)` is executed before (only once).
    auto rate(ChemicalState const& state) const -> real;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Compare two Reaction instances for less than
auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool;

/// Compare two Reaction instances for equality
auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool;

} // namespace Reaktoro
