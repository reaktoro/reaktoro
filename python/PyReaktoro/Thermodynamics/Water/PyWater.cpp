// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PyWater.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateHGK.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzStateWagnerPruss.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {

auto export_WaterConstants() -> void
{
    py::scope().attr("waterMolarMass") = waterMolarMass;
    py::scope().attr("waterCriticalTemperature") = waterCriticalTemperature;
    py::scope().attr("waterCriticalPressure") = waterCriticalPressure;
    py::scope().attr("waterCriticalDensity") = waterCriticalDensity;
    py::scope().attr("waterTriplePointTemperature") = waterTriplePointTemperature;
    py::scope().attr("waterTriplePointPressure") = waterTriplePointPressure;
    py::scope().attr("waterTriplePointDensityLiquid") = waterTriplePointDensityLiquid;
    py::scope().attr("waterTriplePointDensityVapour") = waterTriplePointDensityVapour;
}

auto export_WaterThermoState() -> void
{
    py::class_<WaterThermoState>("WaterThermoState")
        .def_readwrite("temperature", &WaterThermoState::temperature)
        .def_readwrite("volume", &WaterThermoState::volume)
        .def_readwrite("entropy", &WaterThermoState::entropy)
        .def_readwrite("helmholtz", &WaterThermoState::helmholtz)
        .def_readwrite("internal_energy", &WaterThermoState::internal_energy)
        .def_readwrite("enthalpy", &WaterThermoState::enthalpy)
        .def_readwrite("gibbs", &WaterThermoState::gibbs)
        .def_readwrite("cv", &WaterThermoState::cv)
        .def_readwrite("cp", &WaterThermoState::cp)
        .def_readwrite("density", &WaterThermoState::density)
        .def_readwrite("densityT", &WaterThermoState::densityT)
        .def_readwrite("densityP", &WaterThermoState::densityP)
        .def_readwrite("densityTT", &WaterThermoState::densityTT)
        .def_readwrite("densityTP", &WaterThermoState::densityTP)
        .def_readwrite("densityPP", &WaterThermoState::densityPP)
        .def_readwrite("pressure", &WaterThermoState::pressure)
        .def_readwrite("pressureT", &WaterThermoState::pressureT)
        .def_readwrite("pressureD", &WaterThermoState::pressureD)
        .def_readwrite("pressureTT", &WaterThermoState::pressureTT)
        .def_readwrite("pressureTD", &WaterThermoState::pressureTD)
        .def_readwrite("pressureDD", &WaterThermoState::pressureDD)
        ;

    py::def("waterThermoStateHGK", waterThermoStateHGK);
    py::def("waterThermoStateWagnerPruss", waterThermoStateWagnerPruss);
    py::def("waterThermoState", waterThermoState);
}

auto export_WaterHelmholtzState() -> void
{
    py::class_<WaterHelmholtzState>("WaterHelmholtzState")
        .def_readwrite("helmholtz", &WaterHelmholtzState::helmholtz)
        .def_readwrite("helmholtzT", &WaterHelmholtzState::helmholtzT)
        .def_readwrite("helmholtzD", &WaterHelmholtzState::helmholtzD)
        .def_readwrite("helmholtzTT", &WaterHelmholtzState::helmholtzTT)
        .def_readwrite("helmholtzTD", &WaterHelmholtzState::helmholtzTD)
        .def_readwrite("helmholtzDD", &WaterHelmholtzState::helmholtzDD)
        .def_readwrite("helmholtzTTT", &WaterHelmholtzState::helmholtzTTT)
        .def_readwrite("helmholtzTTD", &WaterHelmholtzState::helmholtzTTD)
        .def_readwrite("helmholtzTDD", &WaterHelmholtzState::helmholtzTDD)
        .def_readwrite("helmholtzDDD", &WaterHelmholtzState::helmholtzDDD)
        ;

    py::def("waterHelmholtzStateHGK", waterHelmholtzStateHGK);
    py::def("waterHelmholtzStateWagnerPruss", waterHelmholtzStateWagnerPruss);
}

auto export_WaterElectroState() -> void
{
    py::class_<WaterElectroState>("WaterElectroState")
        .def_readwrite("epsilon", &WaterElectroState::epsilon)
        .def_readwrite("epsilonT", &WaterElectroState::epsilonT)
        .def_readwrite("epsilonP", &WaterElectroState::epsilonP)
        .def_readwrite("epsilonTT", &WaterElectroState::epsilonTT)
        .def_readwrite("epsilonTP", &WaterElectroState::epsilonTP)
        .def_readwrite("epsilonPP", &WaterElectroState::epsilonPP)
        .def_readwrite("bornZ", &WaterElectroState::bornZ)
        .def_readwrite("bornY", &WaterElectroState::bornY)
        .def_readwrite("bornQ", &WaterElectroState::bornQ)
        .def_readwrite("bornN", &WaterElectroState::bornN)
        .def_readwrite("bornU", &WaterElectroState::bornU)
        .def_readwrite("bornX", &WaterElectroState::bornX)
        ;

    py::def("waterElectroStateJohnsonNorton", waterElectroStateJohnsonNorton);
}

auto export_WaterUtils() -> void
{
    py::def("waterDensityHGK", waterDensityHGK);
    py::def("waterDensityWagnerPruss", waterDensityWagnerPruss);
    py::def("waterLiquidDensityHGK", waterLiquidDensityHGK);
    py::def("waterLiquidDensityWagnerPruss", waterLiquidDensityWagnerPruss);
    py::def("waterVaporDensityHGK", waterVaporDensityHGK);
    py::def("waterVaporDensityWagnerPruss", waterVaporDensityWagnerPruss);
    py::def("waterPressureHGK", waterPressureHGK);
    py::def("waterPressureWagnerPruss", waterPressureWagnerPruss);
    py::def("waterSaturatedPressureWagnerPruss", waterSaturatedPressureWagnerPruss);
    py::def("waterSaturatedLiquidDensityWagnerPruss", waterSaturatedLiquidDensityWagnerPruss);
    py::def("waterSaturatedVapourDensityWagnerPruss", waterSaturatedVapourDensityWagnerPruss);
}

auto export_Water() -> void
{
    export_WaterConstants();
    export_WaterThermoState();
    export_WaterHelmholtzState();
    export_WaterElectroState();
    export_WaterUtils();
}

} // namespace Reaktoro

