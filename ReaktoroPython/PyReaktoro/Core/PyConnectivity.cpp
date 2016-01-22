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

#include "PyConnectivity.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Connectivity.hpp>

namespace Reaktoro {

auto export_Connectivity() -> void
{
    py::class_<Connectivity>("Connectivity")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("indicesElementsInSpecies", &Connectivity::indicesElementsInSpecies, py::return_internal_reference<>())
        .def("indicesElementsInPhase", &Connectivity::indicesElementsInPhase, py::return_internal_reference<>())
        .def("indicesSpeciesInPhase", &Connectivity::indicesSpeciesInPhase, py::return_internal_reference<>())
        .def("indicesSpeciesWithElement", &Connectivity::indicesSpeciesWithElement, py::return_internal_reference<>())
        .def("indicesPhasesWithElement", &Connectivity::indicesPhasesWithElement, py::return_internal_reference<>())
        .def("indexPhaseWithSpecies", &Connectivity::indexPhaseWithSpecies)
        ;
}

} // namespace Reaktoro
