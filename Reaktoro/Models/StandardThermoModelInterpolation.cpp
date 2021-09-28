// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "StandardThermoModelInterpolation.hpp"

// Reaktoro includes
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

/// Return a Vec<Param> object containing all Param objects in @p params.
auto extractParams(const StandardThermoModelParamsInterpolation& params) -> Vec<Param>
{
    return {};
}

/// Return a ModelSerializer for given model parameters in @p params.
auto createModelSerializer(const StandardThermoModelParamsInterpolation& params) -> ModelSerializer
{
    return [=]()
    {
        yaml node;
        node["Interpolation"] = params;
        return node;
    };
}

auto StandardThermoModelInterpolation(const StandardThermoModelParamsInterpolation& params) -> StandardThermoModel
{
    const auto& temperatures = params.temperatures;
    const auto& pressures = params.pressures;
    const auto& Pref = params.Pref;

    BilinearInterpolator iG0(temperatures, pressures, params.G0);
    BilinearInterpolator iH0(temperatures, pressures, params.H0);
    BilinearInterpolator iV0(temperatures, pressures, params.V0);
    BilinearInterpolator iVT0(temperatures, pressures, params.VT0);
    BilinearInterpolator iVP0(temperatures, pressures, params.VP0);
    BilinearInterpolator iCp0(temperatures, pressures, params.Cp0);

    auto evalfn = [=](StandardThermoProps& props, real T, real P)
    {
        if(!iG0.empty()) props.G0 = iG0(T, P);
        if(!iH0.empty()) props.H0 = iH0(T, P);
        if(!iV0.empty()) props.V0 = iV0(T, P);
        if(!iVT0.empty()) props.VT0 = iVT0(T, P);
        if(!iVP0.empty()) props.VP0 = iVP0(T, P);
        if(!iCp0.empty()) props.Cp0 = iCp0(T, P);

        props.G0 += props.V0 * (P - Pref);
    };

    return StandardThermoModel(evalfn, extractParams(params), createModelSerializer(params));
}

} // namespace Reaktoro
