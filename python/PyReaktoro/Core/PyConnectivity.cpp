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
