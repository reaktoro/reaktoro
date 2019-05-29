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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Core/Thermo.hpp>

namespace Reaktoro {

void exportThermo(py::module& m)
{
    py::class_<Thermo>(m, "Thermo")
        .def(py::init<const Database&>())
        .def("standardPartialMolarGibbsEnergy", &Thermo::standardPartialMolarGibbsEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarHelmholtzEnergy", &Thermo::standardPartialMolarHelmholtzEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarInternalEnergy", &Thermo::standardPartialMolarInternalEnergy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarEnthalpy", &Thermo::standardPartialMolarEnthalpy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarEntropy", &Thermo::standardPartialMolarEntropy, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarVolume", &Thermo::standardPartialMolarVolume, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarHeatCapacityConstP", &Thermo::standardPartialMolarHeatCapacityConstP, (py::arg("T"), py::arg("P"), "species"))
        .def("standardPartialMolarHeatCapacityConstV", &Thermo::standardPartialMolarHeatCapacityConstV, (py::arg("T"), py::arg("P"), "species"))
        .def("lnEquilibriumConstant", &Thermo::lnEquilibriumConstant, (py::arg("T"), py::arg("P"), "reaction"))
        .def("logEquilibriumConstant", &Thermo::logEquilibriumConstant, (py::arg("T"), py::arg("P"), "reaction"))
        .def("standardPartialMolarGibbsEnergy", &Thermo::standardPartialMolarGibbsEnergy)
        .def("standardPartialMolarHelmholtzEnergy", &Thermo::standardPartialMolarHelmholtzEnergy)
        .def("standardPartialMolarInternalEnergy", &Thermo::standardPartialMolarInternalEnergy)
        .def("standardPartialMolarEnthalpy", &Thermo::standardPartialMolarEnthalpy)
        .def("standardPartialMolarEntropy", &Thermo::standardPartialMolarEntropy)
        .def("standardPartialMolarVolume", &Thermo::standardPartialMolarVolume)
        .def("standardPartialMolarHeatCapacityConstP", &Thermo::standardPartialMolarHeatCapacityConstP)
        .def("standardPartialMolarHeatCapacityConstV", &Thermo::standardPartialMolarHeatCapacityConstV)
        .def("lnEquilibriumConstant", &Thermo::lnEquilibriumConstant)
        .def("logEquilibriumConstant", &Thermo::logEquilibriumConstant)
        ;
}

} // namespace Reaktoro
