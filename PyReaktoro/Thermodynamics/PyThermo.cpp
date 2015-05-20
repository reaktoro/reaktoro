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

#include "PyDatabase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Core/Thermo.hpp>

namespace Reaktoro {

auto export_Thermo() -> void
{
    py::class_<Thermo>("Thermo", py::no_init)
        .def(py::init<const Database&>())
        .def("setTemperatureUnits", &Thermo::setTemperatureUnits)
        .def("setPressureUnits", &Thermo::setPressureUnits)
        .def("standardGibbsEnergy", &Thermo::standardGibbsEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardHelmholtzEnergy", &Thermo::standardHelmholtzEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardInternalEnergy", &Thermo::standardInternalEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardEnthalpy", &Thermo::standardEnthalpy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardEntropy", &Thermo::standardEntropy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardVolume", &Thermo::standardVolume, (py::arg("T"), py::arg("P"), "species"))
        .def("standardHeatCapacity", &Thermo::standardHeatCapacity, (py::arg("T"), py::arg("P"), "species"))
        .def("checkStandardGibbsEnergy", &Thermo::checkStandardGibbsEnergy)
        .def("checkStandardHelmholtzEnergy", &Thermo::checkStandardHelmholtzEnergy)
        .def("checkStandardInternalEnergy", &Thermo::checkStandardInternalEnergy)
        .def("checkStandardEnthalpy", &Thermo::checkStandardEnthalpy)
        .def("checkStandardEntropy", &Thermo::checkStandardEntropy)
        .def("checkStandardVolume", &Thermo::checkStandardVolume)
        .def("checkStandardHeatCapacity", &Thermo::checkStandardHeatCapacity)
        // .def("speciesThermoStateHKF", &Thermo::speciesThermoStateHKF)
        // .def("waterThermoStateHGK", &Thermo::waterThermoStateHGK)
        // .def("waterThermoStateWagnerPruss", &Thermo::waterThermoStateWagnerPruss)
        ;
}

} // namespace Reaktoro
