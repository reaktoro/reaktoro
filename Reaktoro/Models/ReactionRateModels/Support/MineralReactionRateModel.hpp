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
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {

/// The data available for the evaluation of a mineral reaction rate.
struct MineralReactionRateModelArgs
{
    /// The state of the chemical system.
    ChemicalState const& state;

    /// The properties of the chemical system.
    ChemicalProps const& props;

    /// The properties of the aqueous solution.
    AqueousProps const& aprops;

    /// The temperature of the system (in K).
    real const& T;

    /// The pressure of the system (in Pa).
    real const& P;

    /// The pH of the aqueous solution.
    real const& pH;

    /// The saturation index @eq{\Omega = \mathrm{IAP}/K} of the mineral reaction.
    real const& Omega;

    /// The surface area between the mineral and the aqueous solution.
    real const& area;
};

/// The type of functions that calculate rates for mineral dissolution/precipitation reactions.
using MineralReactionRateModel = Model<ReactionRate(MineralReactionRateModelArgs)>;

/// The type of functions that construct a MineralReactionRateModel for a mineral reaction.
/// @param mineral The name of the mineral in the system.
/// @param system The chemical system in which the reactions happen.
using MineralReactionRateModelGenerator = Fn<MineralReactionRateModel(String const& mineral, ChemicalSystem const& system)>;

} // namespace Reaktoro
