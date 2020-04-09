// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

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

void exportWaterConstants(py::module& m)
{
    m.attr("waterMolarMass") = waterMolarMass;
    m.attr("waterCriticalTemperature") = waterCriticalTemperature;
    m.attr("waterCriticalPressure") = waterCriticalPressure;
    m.attr("waterCriticalDensity") = waterCriticalDensity;
    m.attr("waterTriplePointTemperature") = waterTriplePointTemperature;
    m.attr("waterTriplePointPressure") = waterTriplePointPressure;
    m.attr("waterTriplePointDensityLiquid") = waterTriplePointDensityLiquid;
    m.attr("waterTriplePointDensityVapour") = waterTriplePointDensityVapour;
}

void exportWaterThermoState(py::module& m)
{
    py::class_<WaterThermoState>(m, "WaterThermoState")
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

    m.def("waterThermoStateHGK", waterThermoStateHGK);
    m.def("waterThermoStateWagnerPruss", waterThermoStateWagnerPruss);
    m.def("waterThermoState", waterThermoState);
}

void exportWaterHelmholtzState(py::module& m)
{
    py::class_<WaterHelmholtzState>(m, "WaterHelmholtzState")
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

    m.def("waterHelmholtzStateHGK", waterHelmholtzStateHGK);
    m.def("waterHelmholtzStateWagnerPruss", waterHelmholtzStateWagnerPruss);
}

void exportWaterElectroState(py::module& m)
{
    py::class_<WaterElectroState>(m, "WaterElectroState")
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

    m.def("waterElectroStateJohnsonNorton", waterElectroStateJohnsonNorton);
}

void exportWaterUtils(py::module& m)
{
    m.def("waterDensityHGK", waterDensityHGK);
    m.def("waterDensityWagnerPruss", waterDensityWagnerPruss);
    m.def("waterLiquidDensityHGK", waterLiquidDensityHGK);
    m.def("waterLiquidDensityWagnerPruss", waterLiquidDensityWagnerPruss);
    m.def("waterVaporDensityHGK", waterVaporDensityHGK);
    m.def("waterVaporDensityWagnerPruss", waterVaporDensityWagnerPruss);
    m.def("waterPressureHGK", waterPressureHGK);
    m.def("waterPressureWagnerPruss", waterPressureWagnerPruss);
    m.def("waterSaturatedPressureWagnerPruss", waterSaturatedPressureWagnerPruss);
    m.def("waterSaturatedLiquidDensityWagnerPruss", waterSaturatedLiquidDensityWagnerPruss);
    m.def("waterSaturatedVapourDensityWagnerPruss", waterSaturatedVapourDensityWagnerPruss);
}

void exportWater(py::module& m)
{
    exportWaterConstants(m);
    exportWaterThermoState(m);
    exportWaterHelmholtzState(m);
    exportWaterElectroState(m);
    exportWaterUtils(m);
}

} // namespace Reaktoro

