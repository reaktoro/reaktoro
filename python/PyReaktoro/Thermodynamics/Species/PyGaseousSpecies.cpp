// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>

namespace Reaktoro {

void exportGaseousSpecies(py::module& m)
{
    py::class_<GaseousSpecies, Species>(m, "GaseousSpecies")
        .def(py::init<>())
        .def("setCriticalTemperature", &GaseousSpecies::setCriticalTemperature)
        .def("setCriticalPressure", &GaseousSpecies::setCriticalPressure)
        .def("setAcentricFactor", &GaseousSpecies::setAcentricFactor)
        .def("setThermoData", &GaseousSpecies::setThermoData)
        .def("criticalTemperature", &GaseousSpecies::criticalTemperature)
        .def("criticalPressure", &GaseousSpecies::criticalPressure)
        .def("acentricFactor", &GaseousSpecies::acentricFactor)
        .def("thermoData", &GaseousSpecies::thermoData, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
