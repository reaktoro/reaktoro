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
#include <Reaktoro/Core/Reactions.hpp>
#include <Reaktoro/Utils/MineralReactionRateModel.hpp>

namespace Reaktoro {

/// The class used to configure mineral dissolution/precipitation reactions.
class MineralReaction : public GeneralReaction
{
public:
    /// Construct a MineralReaction object with given mineral name.
    /// @param mineral The name of the mineral as found in the database.
    explicit MineralReaction(String const& mineral);

    /// Set the mineral reaction rate model of the reaction.
    auto setRateModel(MineralReactionRateModel const& model) -> MineralReaction&;

    // Consider GeneralReaction::setRateModel methods as well.
    using GeneralReaction::setRateModel;

    /// Set the mineral reaction rate model generator of the reaction.
    /// Use this method to set a MineralReactionRateModelGenerator in case you
    /// need the MineralReactionRateModel of the reaction to be constructed
    /// later, when the chemical system is assembled.
    auto setRateModel(MineralReactionRateModelGenerator const& model_generator) -> MineralReaction&;

    /// Return the name of the mineral.
    auto mineral() const -> String const&;
};

} // namespace Reaktoro
