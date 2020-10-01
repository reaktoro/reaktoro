// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/StandardThermoProps.hpp>
#include <Reaktoro/Core/ReactionThermoProps.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// A class to represent a formation reaction of a chemical species.
/// @ingroup Core
class FormationReaction
{
public:
    /// Construct a default FormationReaction object.
    FormationReaction();

    /// Return a deep copy of this FormationReaction object.
    auto clone() const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new name of the product species in the formation reaction.
    auto withProduct(String product) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reactant species in the formation reaction.
    auto withReactants(Pairs<Species, double> reactants) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new equilibrium constant value.
    /// Note that this method exists for convenience only. Its use results in a
    /// thermodynamic model for this reaction in which its standard enthalpy is
    /// zero. For a more complete thermodynamic model, use method
    /// @ref withReactionThermoPropsFn.
    /// @param lgK0 The equilibrium constant of the reaction (in log base 10)
    auto withEquilibriumConstant(real lgK0) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction thermodynamic model function.
    auto withReactionThermoPropsFn(const ReactionThermoPropsFn& fn) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction thermodynamic model function.
    /// Note that this method exists for convenience only. Consider, for example,
    /// `reaction.with(ReactionThermoModelVantHoff(lgK0, dH0))`, which is equivalent to
    /// `reaction.withReactionThermoPropsFn(ReactionThermoModelVantHoff(lgK0, dH0))`.
    auto with(const ReactionThermoPropsFn& fn) const -> FormationReaction;

    /// Return the name of the product species in the formation reaction.
    auto product() const -> String;

    /// Return the reactant species in the formation reaction.
    auto reactants() const -> const Pairs<Species, double>&;

    /// Return the reaction thermodynamic model function of the formation reaction.
    auto reactionThermoPropsFn() const -> const ReactionThermoPropsFn&;

    /// Return the standard thermodynamic model function of the product species.
    /// This method constructs a standard thermodynamic model function for the
    /// product species using the assigned thermodynamic model of the formation
    /// reaction. An empty model is returned if no reaction thermodynamic model
    /// has been previously assigned.
    auto standardThermoPropsFn() const -> StandardThermoPropsFn;

    /// Return the stoichiometric coefficient of a reactant with given name in the formation reaction.
    auto stoichiometry(String reactant) const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
