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

#include "MineralReactions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {
namespace detail {

/// Convert a vector of MineralReactionRateModel objects to a vector of ReactionRateModel objects.
auto convert(Strings const& minerals, Vec<MineralReactionRateModel> const& models, ChemicalSystem const& system) -> Vec<ReactionRateModel>
{
    errorif(minerals.size() != models.size(), "Expecting same number of minerals and mineral reaction rate models.");

    Vec<ReactionRateModel> ratemodels(models.size());

    auto aprops_ptr = std::make_shared<AqueousProps>(system);

    const auto iaqueousphase = system.phases().indexWithAggregateState(AggregateState::Aqueous);

    for(auto i = 0; i < models.size(); ++i)
    {
        const auto imineral = aprops_ptr->saturationSpecies().indexWithName(minerals[i]);
        const auto imineralphase = system.phases().indexWithName(minerals[i]);
        const auto imineralsurface = system.reactingPhaseInterfaceIndex(imineralphase, imineralphase);

        ratemodels[i] = [=](ChemicalProps const& props) -> ReactionRate
        {
            auto const& aprops = *aprops_ptr;
            auto const& T = props.temperature();
            auto const& P = props.pressure();
            auto const& pH = aprops.pH();
            auto const& Omega = aprops.saturationIndex(imineral);
            auto const& area = props.surfaceArea(imineralsurface);

            const auto args = MineralReactionRateModelArgs{ props, aprops, T, P, pH, Omega, area };

            return models[i](args);
        };
    }

    // Ensure the first rate model updates AqueousProps object before it is used by every mineral rate model!
    ratemodels[0] = [=](ChemicalProps const& props) mutable -> ReactionRate
    {
        aprops_ptr->update(props);
        return ratemodels[0](props);
    };

    return ratemodels;
}

} // namespace detail

MineralReactions::MineralReactions(StringList const& minerals)
: m_minerals(minerals), m_mineral_rate_model_generators(m_minerals.size())
{
}

auto MineralReactions::setRateModel(MineralReactionRateModelGenerator generator) -> void
{
    const auto size = m_mineral_rate_model_generators.size();
    m_mineral_rate_model_generators.assign(size, generator);
}

auto MineralReactions::setRateModel(String mineral, MineralReactionRateModelGenerator generator) -> void
{
    errorif(!generator, "You are trying to specify a non-initialized mineral reaction rate generator.");
    const auto idx = index(m_minerals, mineral);
    errorif(idx >= m_minerals.size(), "You did not specify mineral `", mineral, "` in the list of minerals when creating the MineralReactions object, e.g., `MineralReactions(\"Calcite Dolomite Quartz\")`.");
    m_mineral_rate_model_generators[idx] = generator;
}

auto MineralReactions::operator()(ChemicalSystem const& system) const -> Vec<Reaction>
{
    for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
    {
        const auto imineral = system.species().findWithName(m_minerals[i]);
        errorif(!generator, "You forgot to set a mineral reaction rate model for mineral `", m_minerals[i], "` and maybe for all other minerals as well. Use method MineralReactions::setRateModel to fix this.");
        errorif(imineral >= system.species().size(), "There is no mineral with name `", m_minerals[i], "` in the chemical system.");
    }

    const auto num_minerals = m_minerals.size();

    Vec<MineralReactionRateModel> mineral_reaction_rate_models(num_minerals);
    for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
        mineral_reaction_rate_models[i] = generator(m_minerals[i], system);

    const auto reaction_rate_models = detail::convert(m_minerals, mineral_reaction_rate_models, system);

    Vec<Reaction> reactions(num_minerals);
    for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
    {
        const auto imineral = system.species().index(m_minerals[i]);
        const auto mineralspecies = system.species(imineral);
        reactions[i] = Reaction()
            .withName(m_minerals[i])
            .withEquation({{ {mineralspecies, -1.0} }})
            .withRateModel(reaction_rate_models[i])
            ;
    }

    return reactions;
}

} // namespace Reaktoro
