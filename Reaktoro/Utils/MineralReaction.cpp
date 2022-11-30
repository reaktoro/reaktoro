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

#include "MineralReaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {
namespace detail {

/// Convert a MineralReactionRateModel object to a ReactionRateModel object.
auto convert(String const& mineral, MineralReactionRateModel const& model) -> ReactionRateModel
{
    errorif(!model, "Expecting an initialized MineralReactionRateModel object when converting it to a ReactionRateModel");

    return [=](ChemicalProps const& props) -> ReactionRate
    {
        auto const& aprops = AqueousProps::compute(props);
        auto const& T = props.temperature();
        auto const& P = props.pressure();
        auto const& pH = aprops.pH();
        auto const& Omega = aprops.saturationRatio(mineral);
        auto const& area = props.surfaceArea(mineral);

        const auto args = MineralReactionRateModelArgs{ props, aprops, T, P, pH, Omega, area };

        // Evaluate the mineral reaction rate model and switch sign because
        // Reaktoro's convention for reaction rate is positive when the
        // reaction proceeds from left to right (and this is how the mineral
        // reaction is represented, with the mineral on the left side, its
        // dissolution rate, from left to right, should be positive).
        return -model(args);
    };
}

} // namespace detail

MineralReaction::MineralReaction(String const& mineral)
: GeneralReaction(mineral)
{}

auto MineralReaction::setRateModel(MineralReactionRateModel const& model) -> MineralReaction&
{
    ReactionRateModel converted = detail::convert(mineral(), model);
    GeneralReaction::setRateModel(converted);
    return *this;
}

auto MineralReaction::setRateModel(MineralReactionRateModelGenerator const& model_generator) -> MineralReaction&
{
    const auto mineralname = mineral();
    ReactionRateModelGenerator converted = [=](ReactionRateModelGeneratorArgs args)
    {
        MineralReactionRateModel model = model_generator(args);
        return detail::convert(mineralname, model);
    };
    GeneralReaction::setRateModel(converted);
    return *this;
}

auto MineralReaction::mineral() const -> String const&
{
    return GeneralReaction::name();
}

} // namespace Reaktoro
