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

// Get a persistent AqueousProps object that corresponds to a given ChemicalSystem object.
auto getAqueousProps(ChemicalSystem const& system) -> AqueousProps&
{
    const auto id = system.id();
    thread_local Map<Index, AqueousProps> aprops;
    if(auto it = aprops.find(id); it != aprops.end())
        return it->second;
    const auto [it, _] = aprops.emplace(id, AqueousProps(system));
    return it->second;
}

/// Convert a vector of MineralReactionRateModel objects to a vector of ReactionRateModel objects.
auto convert(Strings const& minerals, Vec<MineralReactionRateModel> const& models) -> Vec<ReactionRateModel>
{
    errorif(minerals.size() != models.size(), "Expecting same number of minerals and mineral reaction rate models.");

    Vec<ReactionRateModel> ratemodels(models.size());

    for(auto i = 0; i < models.size(); ++i)
    {
        ratemodels[i] = [=](ChemicalProps const& props) -> ReactionRate
        {
            auto const& system = props.system();
            auto const& aprops = getAqueousProps(system); // get the AqueousProps object corresponding to the ChemicalSystem object in `props`, which will be reused by all mineral reactions!

            thread_local auto const imineral = aprops.saturationSpecies().indexWithName(minerals[i]);
            thread_local auto const imineralphase = system.phases().indexWithName(minerals[i]);
            thread_local auto const imineralsurface = system.surfaces().indexWithPhases(imineralphase, imineralphase);

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
        getAqueousProps(props.system()).update(props); // update the AqueousProps object corresponding to the ChemicalSystem object in `props`; note this is done once for the first mineral rate model, and updated state reused by all other minerals for performance reasons!
        return ratemodels[0](props);
    };

    return ratemodels;
}

} // namespace detail

MineralReactions::MineralReactions(StringList const& minerals)
: m_minerals(minerals), m_mineral_rate_model_generators(m_minerals.size())
{
}

auto MineralReactions::setRateModel(MineralReactionRateModelGenerator const& generator) -> void
{
    const auto size = m_mineral_rate_model_generators.size();
    m_mineral_rate_model_generators.assign(size, generator);
}

auto MineralReactions::setRateModel(String const& mineral, MineralReactionRateModelGenerator const& generator) -> void
{
    errorif(!generator, "You are trying to specify a non-initialized mineral reaction rate generator.");
    const auto idx = index(m_minerals, mineral);
    errorif(idx >= m_minerals.size(), "You did not specify mineral `", mineral, "` in the list of minerals when creating the MineralReactions object, e.g., `MineralReactions(\"Calcite Dolomite Quartz\")`.");
    m_mineral_rate_model_generators[idx] = generator;
}

auto MineralReactions::operator()(PhaseList const& phases) const -> Vec<Reaction>
{
    auto const& species = phases.species();

    for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
    {
        const auto imineral = species.findWithName(m_minerals[i]);
        errorif(!generator, "You forgot to set a mineral reaction rate model for mineral `", m_minerals[i], "` and maybe for all other minerals as well. Use method MineralReactions::setRateModel to fix this.");
        errorif(imineral >= species.size(), "There is no mineral with name `", m_minerals[i], "` in the chemical system.");
    }

    const auto num_minerals = m_minerals.size();

    Vec<MineralReactionRateModel> mineral_reaction_rate_models(num_minerals);
    for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
        mineral_reaction_rate_models[i] = generator(m_minerals[i], phases);

    const auto reaction_rate_models = detail::convert(m_minerals, mineral_reaction_rate_models);

    Vec<Reaction> reactions(num_minerals);
    for(auto i = 0; i < num_minerals; ++i)
    {
        const auto imineral = species.index(m_minerals[i]);
        const auto mineralspecies = species[imineral];
        reactions[i] = Reaction()
            .withName(m_minerals[i])
            .withEquation({{ {mineralspecies, -1.0} }})
            .withRateModel(reaction_rate_models[i])
            ;
    }

    return reactions;
}

} // namespace Reaktoro
