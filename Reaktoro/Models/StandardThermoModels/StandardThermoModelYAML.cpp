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

#include "StandardThermoModelYAML.hpp"

// Reaktoro includes
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

auto StandardThermoModelYAML(const yaml& node) -> StandardThermoModel
{
    errorif(!(node.IsMap() && node.size() == 1),
        "Expecting the following yaml node to contain a single "
        "StandardThermoModel child:\n", node.repr());

    auto model = node.begin()->first.as<String>();
    yaml params = node.begin()->second;

    if(model == "Constant")
        return StandardThermoModelConstant(params);
    if(model == "HKF")
        return StandardThermoModelHKF(params);
    if(model == "HollandPowell")
        return StandardThermoModelHollandPowell(params);
    if(model == "Interpolation")
        return StandardThermoModelInterpolation(params);
    if(model == "MaierKelley")
        return StandardThermoModelMaierKelley(params);
    if(model == "MineralHKF")
        return StandardThermoModelMineralHKF(params);
    if(model == "WaterHKF")
        return StandardThermoModelWaterHKF(params);
    if(model == "Nasa")
        return StandardThermoModelNasa(params);
    errorif(true, "Cannot create a StandardThermoModel with "
        "unsupported model name `", model, "` in yaml node:\n", node.repr());

    return {};
}

} // namespace Reaktoro
