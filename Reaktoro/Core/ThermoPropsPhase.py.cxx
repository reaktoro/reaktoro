// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ThermoPropsPhase.hpp>
using namespace Reaktoro;

void exportThermoPropsPhase(py::module& m)
{
    auto update1 = [](ThermoPropsPhase& self, const real& T, const real& P)
    {
        self.update(T, P);
    };

    auto update2 = [](ThermoPropsPhase& self, const real& T, const real& P, Wrt<real&> wrtvar)
    {
        self.update(T, P, wrtvar);
    };

    py::class_<ThermoPropsPhase>(m, "ThermoPropsPhase")
        .def(py::init<const Phase&>())
        .def("update", update1)
        .def("update", update2)
        .def("phase", &ThermoPropsPhase::phase, py::return_value_policy::reference_internal)
        .def("data", &ThermoPropsPhase::data, py::return_value_policy::reference_internal)
        .def("temperature", &ThermoPropsPhase::temperature)
        .def("pressure", &ThermoPropsPhase::pressure)
        .def("standardVolumes", &ThermoPropsPhase::standardVolumes, py::return_value_policy::reference_internal)
        .def("standardGibbsEnergies", &ThermoPropsPhase::standardGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardEnthalpies", &ThermoPropsPhase::standardEnthalpies, py::return_value_policy::reference_internal)
        .def("standardEntropies", &ThermoPropsPhase::standardEntropies)
        .def("standardInternalEnergies", &ThermoPropsPhase::standardInternalEnergies)
        .def("standardHelmholtzEnergies", &ThermoPropsPhase::standardHelmholtzEnergies)
        .def("standardHeatCapacitiesConstP", &ThermoPropsPhase::standardHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardHeatCapacitiesConstV", &ThermoPropsPhase::standardHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        ;
}
