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
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Models/ReactionRateModels/Support/MineralReactionRateModel.hpp>

namespace Reaktoro {

/// The class used to configure mineral dissolution/precipitation reactions.
class MineralReactions
{
public:
    /// Construct a MineralReactions object with given mineral names.
    MineralReactions(StringList const& minerals);

    /// Set a common mineral reaction rate model generator for all minerals.
    auto setRateModel(MineralReactionRateModelGenerator const& generator) -> MineralReactions&;

    /// Set a mineral reaction rate model generator for a specific mineral.
    auto setRateModel(String const& mineral, MineralReactionRateModelGenerator const& generator) -> MineralReactions&;

    /// Convert this MineralReactions object into a vector of Reaction objects (for compatibility with ReactionGenerator).
    /// @param species The species composing the chemical system where the reactions will take place.
    auto operator()(SpeciesList const& species) const -> Vec<Reaction>;

private:
    /// The names of the minerals
    Strings m_minerals;

    /// The reaction rate model generators for each mineral reaction
    Vec<MineralReactionRateModelGenerator> m_mineral_rate_model_generators;
};

} // namespace Reaktoro
