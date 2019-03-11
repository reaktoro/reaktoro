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
