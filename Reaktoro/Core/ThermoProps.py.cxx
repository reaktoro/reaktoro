// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ThermoProps.hpp>
using namespace Reaktoro;

void exportThermoProps(py::module& m)
{
    auto update1 = [](ThermoProps& self, const real& T, const real& P)
    {
        self.update(T, P);
    };

    auto update2 = [](ThermoProps& self, const real& T, const real& P, Wrt<real&> wrtvar)
    {
        self.update(T, P, wrtvar);
    };

    py::class_<ThermoProps>(m, "ThermoProps")
        .def(py::init<const ChemicalSystem&>())
        .def("update", update1)
        .def("update", update2)
        .def("system", &ThermoProps::system, py::return_value_policy::reference_internal)
        .def("phaseProps", &ThermoProps::phaseProps, py::return_value_policy::reference_internal)
        .def("temperature", &ThermoProps::temperature, py::return_value_policy::reference_internal)
        .def("pressure", &ThermoProps::pressure, py::return_value_policy::reference_internal)
        .def("standardVolumes", &ThermoProps::standardVolumes, py::return_value_policy::reference_internal)
        .def("standardGibbsEnergies", &ThermoProps::standardGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardEnthalpies", &ThermoProps::standardEnthalpies, py::return_value_policy::reference_internal)
        .def("standardEntropies", &ThermoProps::standardEntropies)
        .def("standardInternalEnergies", &ThermoProps::standardInternalEnergies)
        .def("standardHelmholtzEnergies", &ThermoProps::standardHelmholtzEnergies)
        .def("standardHeatCapacitiesConstP", &ThermoProps::standardHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardHeatCapacitiesConstV", &ThermoProps::standardHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        ;
}
