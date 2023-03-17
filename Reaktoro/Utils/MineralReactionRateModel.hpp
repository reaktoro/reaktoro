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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/ReactionRate.hpp>
#include <Reaktoro/Core/ReactionRateModel.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {

/// The data provided to a MineralReactionRateModel to evaluate the rate of a mineral reaction.
/// @see MineralReactionRateModel
struct MineralReactionRateModelArgs
{
    /// The properties of the chemical system.
    ChemicalProps const& props;

    /// The properties of the aqueous solution.
    AqueousProps const& aprops;

    /// The temperature of the chemical system (in K).
    real const& T;

    /// The pressure of the chemical system (in Pa).
    real const& P;

    /// The pH of the aqueous solution.
    real const& pH;

    /// The saturation ratio @eq{\Omega = \mathrm{IAP}/K} of the mineral reaction.
    real const& Omega;

    /// The surface area between the mineral and the aqueous solution.
    real const& area;
};

/// The type of functions that calculate rates for mineral dissolution/precipitation reactions.
/// The sign convention for a mineral reaction rate is negative if dissolving, positive if precipitating.
/// @param args The data provided to evaluate the mineral reaction rate.
using MineralReactionRateModel = Model<ReactionRate(MineralReactionRateModelArgs args)>;

/// The type of functions that construct a MineralReactionRateModel for a mineral reaction.
/// @param args The data provided to construct a MineralReactionRateModel object.
/// @see MineralReaction
using MineralReactionRateModelGenerator = Fn<MineralReactionRateModel(ReactionRateModelGeneratorArgs args)>;

} // namespace Reaktoro
