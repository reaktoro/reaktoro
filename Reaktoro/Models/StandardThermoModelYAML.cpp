// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Models/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Serialization/Serialization.yaml.hpp>

namespace Reaktoro {

auto StandardThermoModelYAML(const String& model, const yaml& params) -> StandardThermoModel
{
    if(model == "HKF")
        return StandardThermoModelHKF(params);
    if(model == "HollandPowell")
        return StandardThermoModelHollandPowell(params);
    if(model == "MaierKelley")
        return StandardThermoModelMaierKelley(params);
    if(model == "MineralHKF")
        return StandardThermoModelMineralHKF(params);
    if(model == "WaterHKF")
        return StandardThermoModelWaterHKF(params);
    errorif(true, "Cannot create StandardThermoModel with not supported model name `", model, "`.");
}

} // namespace Reaktoro
