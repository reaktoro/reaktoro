// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "ActivityModel.hpp"

namespace Reaktoro {

auto chain(Vec<ActivityModelGenerator> const& models) -> ActivityModelGenerator
{
    ActivityModelGenerator chained_model = [=](SpeciesList const& species)
    {
        const Vec<ActivityModel> activity_models = vectorize(models, RKT_LAMBDA(model, model(species)));

        ActivityModel chained_activity_model = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            for(const auto& fn : activity_models)
                fn(props, args);
        };

        return chained_activity_model;
    };

    return chained_model;
}

auto chain(ActivityModelGenerator const& model) -> ActivityModelGenerator
{
    return model;
}

} // namespace Reaktoro
