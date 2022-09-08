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
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/Rate.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Reaction;

/// The type of functions for calculation of reaction rates (in mol/s).
/// @param state The state of the chemical system
/// @return The rate of the reaction (in mol/s)
/// @see Reaction
/// @ingroup Core
using ReactionRateModel = Model<Rate(ChemicalState const& state)>;

/// The type of functions that construct a ReactionRateModel for a reaction.
/// @param reaction The reaction for which the reaction rate model is constructed.
/// @param system The chemical system in which the reactions happen.
/// @ingroup Core
using ReactionRateModelGenerator = Fn<ReactionRateModel(Reaction const& reaction, ChemicalSystem const& system)>;

} // namespace Reaktoro
