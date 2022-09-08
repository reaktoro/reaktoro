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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelConstLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/StandardThermoModels/ReactionStandardThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelNasa.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Models/StandardThermoModels/StandardVolumeModelConstant.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>
using namespace Reaktoro;

void exportSerializationModelsYAML(py::module& m)
{
    //======================================================================
    // StandardThermoModelParams Types
    //======================================================================
    py::implicitly_convertible<yaml, StandardThermoModelParamsConstant>();
    py::implicitly_convertible<StandardThermoModelParamsConstant, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsHKF>();
    py::implicitly_convertible<StandardThermoModelParamsHKF, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsHollandPowell>();
    py::implicitly_convertible<StandardThermoModelParamsHollandPowell, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsInterpolation>();
    py::implicitly_convertible<StandardThermoModelParamsInterpolation, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsMaierKelley>();
    py::implicitly_convertible<StandardThermoModelParamsMaierKelley, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsMineralHKF>();
    py::implicitly_convertible<StandardThermoModelParamsMineralHKF, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsNasa>();
    py::implicitly_convertible<StandardThermoModelParamsNasa, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModelParamsWaterHKF>();
    py::implicitly_convertible<StandardThermoModelParamsWaterHKF, yaml>();

    //======================================================================
    // ReactionStandardThermoModelParams Types
    //======================================================================
    py::implicitly_convertible<yaml, ReactionStandardThermoModelParamsConstLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsConstLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionStandardThermoModelParamsGemsLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsGemsLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionStandardThermoModelParamsPhreeqcLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsPhreeqcLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionStandardThermoModelParamsVantHoff>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsVantHoff, yaml>();

    //======================================================================
    // StandardVolumeModelParams Types
    //======================================================================
    py::implicitly_convertible<yaml, StandardVolumeModelParamsConstant>();
    py::implicitly_convertible<StandardVolumeModelParamsConstant, yaml>();

    //======================================================================
    // ReactionRateModelParams Types
    //======================================================================
    py::implicitly_convertible<yaml, ReactionRateModelParamsPalandriKharaka>();
    py::implicitly_convertible<ReactionRateModelParamsPalandriKharaka, yaml>();
}
