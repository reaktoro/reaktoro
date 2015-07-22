// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "PyThermoProperties.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ThermoProperties.hpp>

namespace Reaktoro {

auto export_ThermoProperties() -> void
{
    py::class_<ThermoProperties>("ThermoProperties")
        .def(py::init<>())
        .def(py::init<unsigned>())
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

    py::class_<SpeciesThermoProperties>("SpeciesThermoProperties")
        .def(py::init<>())
        .def("temperature", &SpeciesThermoProperties::temperature)
        .def("pressure", &SpeciesThermoProperties::pressure)
        .def("standardPartialMolarGibbsEnergy", &SpeciesThermoProperties::standardPartialMolarGibbsEnergy)
        .def("standardPartialMolarEnthalpy", &SpeciesThermoProperties::standardPartialMolarEnthalpy)
        .def("standardPartialMolarVolume", &SpeciesThermoProperties::standardPartialMolarVolume)
        .def("standardPartialMolarEntropy", &SpeciesThermoProperties::standardPartialMolarEntropy)
        .def("standardPartialMolarInternalEnergy", &SpeciesThermoProperties::standardPartialMolarInternalEnergy)
        .def("standardPartialMolarHelmholtzEnergy", &SpeciesThermoProperties::standardPartialMolarHelmholtzEnergy)
        .def("standardPartialMolarHeatCapacityConstP", &SpeciesThermoProperties::standardPartialMolarHeatCapacityConstP)
        .def("standardPartialMolarHeatCapacityConstV", &SpeciesThermoProperties::standardPartialMolarHeatCapacityConstV)
        ;
}

} // namespace Reaktoro
