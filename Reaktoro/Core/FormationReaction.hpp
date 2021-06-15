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

    /// Return a duplicate of this FormationReaction object with new reactant species in the formation reaction.
    auto withReactants(Pairs<Species, double> reactants) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction
    /// thermodynamic model function.
    ///
    /// This method sets a constant equilibrium constant for the formation reaction.
    /// It also sets the standard molar enthalpy and volume of the reaction to zero.
    ///
    /// For a more advanced thermodynamic model setup, use methods
    /// @ref withProductStandardVolumeModel and @ref withReactionThermoModel.
    ///
    /// @param lgK0 The equilibrium constant of the reaction (in log base 10)
    auto withEquilibriumConstant(Param lgK0) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new standard
    /// molar volume function for the product species.
    ///
    /// The standard molar volume of the product species is needed to enable
    /// the calculation of the standard molar volume of the reaction.
    ///
    /// @warning Ensure the standard molar volume of the product species is set
    /// using either @ref withProductStandardVolume or @ref withProductStandardVolumeModel.
    /// Calling @ref withEquilibriumConstant also sets the standard molar volume of the
    /// product species, but to zero.
    ///
    /// @param V0p The constant standard molar volume of the product species (in m3/mol).
    auto withProductStandardVolume(Param V0p) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new standard
    /// molar volume function for the product species.
    ///
    /// The standard molar volume of the product species is needed to enable
    /// the calculation of the standard molar volume of the reaction.
    ///
    /// @warning Ensure the standard molar volume of the product species is set
    /// using either @ref withProductStandardVolume or @ref withProductStandardVolumeModel.
    /// Calling @ref withEquilibriumConstant also sets the standard molar volume of the
    /// product species, but to zero.
    ///
    /// @param fn The standard molar volume model of the product species (in m3/mol).
    auto withProductStandardVolumeModel(Model<real(real,real)> fn) const -> FormationReaction;

    /// Return a duplicate of this FormationReaction object with new reaction thermodynamic model function.
    auto withReactionThermoModel(const ReactionThermoModel& fn) const -> FormationReaction;

    /// Return the reactant species in the formation reaction.
    auto reactants() const -> const Pairs<Species, double>&;

    /// Return the stoichiometric coefficient of a reactant with given name in the formation reaction.
    auto stoichiometry(String reactant) const -> double;

    /// Return the standard molar volume function of the product species.
    auto productStandardVolumeModel() const -> const Model<real(real,real)>&;

    /// Return the reaction thermodynamic model function of the formation reaction.
    auto reactionThermoModel() const -> const ReactionThermoModel&;

    /// Construct the standard thermodynamic model function of the product species.
    ///
    /// This method constructs a standard thermodynamic model function for the
    /// product species using the assigned thermodynamic model of the formation
    /// reaction. An empty model is returned if no reaction thermodynamic model
    /// has been previously assigned.
    ///
    /// @warning This method will throw a runtime error if methods
    /// @ref withProductStandardVolumeModel and @ref withReactionThermoModel
    /// have not been invoked.
    auto createStandardThermoModel() const -> StandardThermoModel;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
