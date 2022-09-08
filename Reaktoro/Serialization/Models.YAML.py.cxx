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
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModels/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>
using namespace Reaktoro;

void exportSerializationModelsYAML(py::module& m)
{
    py::implicitly_convertible<yaml, ReactionThermoModelParamsConstLgK>();
    py::implicitly_convertible<ReactionThermoModelParamsConstLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionThermoModelParamsGemsLgK>();
    py::implicitly_convertible<ReactionThermoModelParamsGemsLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionThermoModelParamsPhreeqcLgK>();
    py::implicitly_convertible<ReactionThermoModelParamsPhreeqcLgK, yaml>();

    py::implicitly_convertible<yaml, ReactionThermoModelParamsVantHoff>();
    py::implicitly_convertible<ReactionThermoModelParamsVantHoff, yaml>();

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

    py::implicitly_convertible<yaml, StandardThermoModelParamsWaterHKF>();
    py::implicitly_convertible<StandardThermoModelParamsWaterHKF, yaml>();
}
