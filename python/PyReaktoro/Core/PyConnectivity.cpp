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
#include <Reaktoro/Core/Connectivity.hpp>

namespace Reaktoro {

void exportConnectivity(py::module& m)
{
    py::class_<Connectivity>(m, "Connectivity")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("indicesElementsInSpecies", &Connectivity::indicesElementsInSpecies, py::return_value_policy::reference_internal)
        .def("indicesElementsInPhase", &Connectivity::indicesElementsInPhase, py::return_value_policy::reference_internal)
        .def("indicesSpeciesInPhase", &Connectivity::indicesSpeciesInPhase, py::return_value_policy::reference_internal)
        .def("indicesSpeciesWithElement", &Connectivity::indicesSpeciesWithElement, py::return_value_policy::reference_internal)
        .def("indicesPhasesWithElement", &Connectivity::indicesPhasesWithElement, py::return_value_policy::reference_internal)
        .def("indexPhaseWithSpecies", &Connectivity::indexPhaseWithSpecies)
        ;
}

} // namespace Reaktoro
