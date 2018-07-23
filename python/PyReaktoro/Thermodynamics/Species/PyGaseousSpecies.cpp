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

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>

// PyReator includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_GaseousSpecies() -> void
{
    py::class_<GaseousSpecies, py::bases<Species>>("GaseousSpecies")
        .def(py::init<>())
        .def("setCriticalTemperature", &GaseousSpecies::setCriticalTemperature)
        .def("setCriticalPressure", &GaseousSpecies::setCriticalPressure)
        .def("setAcentricFactor", &GaseousSpecies::setAcentricFactor)
        .def("setThermoData", &GaseousSpecies::setThermoData)
        .def("criticalTemperature", &GaseousSpecies::criticalTemperature)
        .def("criticalPressure", &GaseousSpecies::criticalPressure)
        .def("acentricFactor", &GaseousSpecies::acentricFactor)
        .def("thermoData", &GaseousSpecies::thermoData, py::return_internal_reference<>())
        ;

    export_std_vector<GaseousSpecies>("GaseousSpeciesVector");
}

} // namespace Reaktoro
