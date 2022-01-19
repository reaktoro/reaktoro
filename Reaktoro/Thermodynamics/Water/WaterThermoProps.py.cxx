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
#include <Reaktoro/Thermodynamics/Water/WaterThermoProps.hpp>
using namespace Reaktoro;

void exportWaterThermoProps(py::module& m)
{
    py::class_<WaterThermoProps>(m, "WaterThermoProps")
        .def(py::init<>())
        .def_readwrite("temperature", &WaterThermoProps::temperature)
        .def_readwrite("volume", &WaterThermoProps::volume)
        .def_readwrite("entropy", &WaterThermoProps::entropy)
        .def_readwrite("helmholtz", &WaterThermoProps::helmholtz)
        .def_readwrite("internal_energy", &WaterThermoProps::internal_energy)
        .def_readwrite("enthalpy", &WaterThermoProps::enthalpy)
        .def_readwrite("gibbs", &WaterThermoProps::gibbs)
        .def_readwrite("cv", &WaterThermoProps::cv)
        .def_readwrite("cp", &WaterThermoProps::cp)
        .def_readwrite("density", &WaterThermoProps::density)
        .def_readwrite("densityT", &WaterThermoProps::densityT)
        .def_readwrite("densityP", &WaterThermoProps::densityP)
        .def_readwrite("densityTT", &WaterThermoProps::densityTT)
        .def_readwrite("densityTP", &WaterThermoProps::densityTP)
        .def_readwrite("densityPP", &WaterThermoProps::densityPP)
        .def_readwrite("pressure", &WaterThermoProps::pressure)
        .def_readwrite("pressureT", &WaterThermoProps::pressureT)
        .def_readwrite("pressureD", &WaterThermoProps::pressureD)
        .def_readwrite("pressureTT", &WaterThermoProps::pressureTT)
        .def_readwrite("pressureDD", &WaterThermoProps::pressureDD)
        .def_readwrite("pressureTD", &WaterThermoProps::pressureTD)
        ;
}
