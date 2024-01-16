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

#include "StandardThermoModelFromData.hpp"

// Reaktoro includes
#include <Reaktoro/Models/StandardThermoModels.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>

namespace Reaktoro {

auto StandardThermoModelFromData(Data const& data) -> StandardThermoModel
{
    errorif(!data.isDict(),
        "Expecting a dictionary containing a single key-value pair in the Data object:\n", data.repr());

    const auto model = data.asDict();

    errorif(model.size() != 1,
        "Expecting only one key-value pair in the Data object:\n", data.repr());

    auto const& name = model.front().first;
    auto const& params = model.front().second;

    if(name == "Constant")
        return StandardThermoModelConstant(params.as<StandardThermoModelParamsConstant>());
    if(name == "ExtendedUNIQUAC")
        return StandardThermoModelExtendedUNIQUAC(params.as<StandardThermoModelParamsExtendedUNIQUAC>());
    if(name == "HKF")
        return StandardThermoModelHKF(params.as<StandardThermoModelParamsHKF>());
    if(name == "HollandPowell")
        return StandardThermoModelHollandPowell(params.as<StandardThermoModelParamsHollandPowell>());
    if(name == "Interpolation")
        return StandardThermoModelInterpolation(params.as<StandardThermoModelParamsInterpolation>());
    if(name == "MaierKelley")
        return StandardThermoModelMaierKelley(params.as<StandardThermoModelParamsMaierKelley>());
    if(name == "MineralHKF")
        return StandardThermoModelMineralHKF(params.as<StandardThermoModelParamsMineralHKF>());
    if(name == "WaterHKF")
        return StandardThermoModelWaterHKF(params.as<StandardThermoModelParamsWaterHKF>());
    if(name == "Nasa")
        return StandardThermoModelNasa(params.as<StandardThermoModelParamsNasa>());

    errorif(true, "Cannot create a StandardThermoModel object with "
        "unsupported model name `", name, "` in Data object:\n", data.repr());

    return {};
}

} // namespace Reaktoro
