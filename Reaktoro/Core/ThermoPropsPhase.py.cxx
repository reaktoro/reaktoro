// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
        .def("phase", &ThermoPropsPhase::phase, return_internal_ref)
        .def("data", &ThermoPropsPhase::data, return_internal_ref)
        .def("temperature", &ThermoPropsPhase::temperature)
        .def("pressure", &ThermoPropsPhase::pressure)
        .def("speciesStandardVolumes", &ThermoPropsPhase::speciesStandardVolumes, return_internal_ref)
        .def("speciesStandardVolumesT", &ThermoPropsPhase::speciesStandardVolumesT, return_internal_ref)
        .def("speciesStandardVolumesP", &ThermoPropsPhase::speciesStandardVolumesP, return_internal_ref)
        .def("speciesStandardGibbsEnergies", &ThermoPropsPhase::speciesStandardGibbsEnergies, return_internal_ref)
        .def("speciesStandardEnthalpies", &ThermoPropsPhase::speciesStandardEnthalpies, return_internal_ref)
        .def("speciesStandardEntropies", &ThermoPropsPhase::speciesStandardEntropies)
        .def("speciesStandardInternalEnergies", &ThermoPropsPhase::speciesStandardInternalEnergies)
        .def("speciesStandardHelmholtzEnergies", &ThermoPropsPhase::speciesStandardHelmholtzEnergies)
        .def("speciesStandardHeatCapacitiesConstP", &ThermoPropsPhase::speciesStandardHeatCapacitiesConstP, return_internal_ref)
        .def("speciesStandardHeatCapacitiesConstV", &ThermoPropsPhase::speciesStandardHeatCapacitiesConstV, return_internal_ref)
        ;
}
