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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>

namespace Reaktoro {

void exportThermoProperties(py::module& m)
{
    py::class_<ThermoProperties>(m, "ThermoProperties")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("update", &ThermoProperties::update)
        .def("temperature", &ThermoProperties::temperature)
        .def("pressure", &ThermoProperties::pressure)
        .def("standardPartialMolarGibbsEnergies", &ThermoProperties::standardPartialMolarGibbsEnergies)
        .def("standardPartialMolarEnthalpies", &ThermoProperties::standardPartialMolarEnthalpies)
        .def("standardPartialMolarVolumes", &ThermoProperties::standardPartialMolarVolumes)
        .def("standardPartialMolarEntropies", &ThermoProperties::standardPartialMolarEntropies)
        .def("standardPartialMolarInternalEnergies", &ThermoProperties::standardPartialMolarInternalEnergies)
        .def("standardPartialMolarHelmholtzEnergies", &ThermoProperties::standardPartialMolarHelmholtzEnergies)
        .def("standardPartialMolarHeatCapacitiesConstP", &ThermoProperties::standardPartialMolarHeatCapacitiesConstP)
        .def("standardPartialMolarHeatCapacitiesConstV", &ThermoProperties::standardPartialMolarHeatCapacitiesConstV)
        ;
}

} // namespace Reaktoro
