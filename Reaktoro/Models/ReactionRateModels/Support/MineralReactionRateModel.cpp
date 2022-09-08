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

#include "MineralReactionRateModel.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {

auto convert(Strings const& minerals, Vec<MineralReactionRateModel> const& models, Database const& database, PhaseList& phases) -> Vec<ReactionRateModel>
{
    errorif(minerals.size() != models.size(), "Expecting same number of minerals and mineral reaction rate models.");

    Vec<ReactionRateModel> ratemodels(models.size());

    ChemicalSystem system(database, phases);

    auto aprops_ptr = std::make_shared<AqueousProps>(system);

    const auto iaqueousphase = phases.indexWithAggregateState(AggregateState::Aqueous);

    for(auto i = 0; i < models.size(); ++i)
    {
        const auto imineral = aprops_ptr->saturationSpecies().indexWithName(minerals[i]);
        const auto imineralphase = system.phases().indexWithName(minerals[i]);

        ratemodels[i] = [=](ChemicalState const& state)
        {
            auto const& props = state.props();
            auto const& aprops = *aprops_ptr;
            auto const& T = props.temperature();
            auto const& P = props.pressure();
            auto const& pH = aprops.pH();
            auto const& Omega = aprops.saturationIndex(imineral);
            auto const& area = state.surfaceArea(iaqueousphase, imineralphase);

            const auto args = MineralReactionRateArgs{ state, props, aprops, T, P, pH, Omega, area };

            return models[i](args);
        };
    }

    // Ensure the first rate model updates AqueousProps object before it is used by every mineral rate model!
    ratemodels[0] = [=](ChemicalState const& state)
    {
        aprops_ptr->update(state);
        return ratemodels[0](state);
    };

    return ratemodels;
}

} // namespace Reaktoro
