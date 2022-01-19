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

#include "StandardVolumeModelConstant.hpp"

// Reaktoro includes
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(const StandardVolumeModelParamsConstant& params) -> Vec<Param>
{
    const auto& [V0] = params;
    return {V0};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardVolumeModelParamsConstant& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["Constant"] = params;
        return node;
    };
}

auto StandardVolumeModelConstant(const StandardVolumeModelParamsConstant& params) -> StandardVolumeModel
{
    auto calcfn = [=](real T, real P)
    {
        return params.V0;
    };

    return StandardVolumeModel(calcfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
