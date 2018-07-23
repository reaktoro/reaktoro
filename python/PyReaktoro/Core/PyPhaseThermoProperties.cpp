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

#include "PyPhaseThermoProperties.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseThermoProperties.hpp>

namespace Reaktoro {

auto export_PhaseThermoProperties() -> void
{
    py::class_<PhaseThermoProperties>("PhaseThermoProperties")
        .def(py::init<>())
        .def(py::init<const Phase&>())
        .def("update", &PhaseThermoProperties::update)
        .def("temperature", &PhaseThermoProperties::temperature)
        .def("pressure", &PhaseThermoProperties::pressure)
        .def("standardPartialMolarGibbsEnergies", &PhaseThermoProperties::standardPartialMolarGibbsEnergies)
        .def("standardPartialMolarEnthalpies", &PhaseThermoProperties::standardPartialMolarEnthalpies)
        .def("standardPartialMolarVolumes", &PhaseThermoProperties::standardPartialMolarVolumes)
        .def("standardPartialMolarEntropies", &PhaseThermoProperties::standardPartialMolarEntropies)
        .def("standardPartialMolarInternalEnergies", &PhaseThermoProperties::standardPartialMolarInternalEnergies)
        .def("standardPartialMolarHelmholtzEnergies", &PhaseThermoProperties::standardPartialMolarHelmholtzEnergies)
        .def("standardPartialMolarHeatCapacitiesConstP", &PhaseThermoProperties::standardPartialMolarHeatCapacitiesConstP)
        .def("standardPartialMolarHeatCapacitiesConstV", &PhaseThermoProperties::standardPartialMolarHeatCapacitiesConstV)
        ;
}

} // namespace Reaktoro
