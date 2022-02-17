// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/PhaseList.hpp>
using namespace Reaktoro;

void exportPhaseList(py::module& m)
{
    py::class_<PhaseList>(m, "PhaseList")
        .def(py::init<>())
        .def(py::init<const Vec<Phase>&>())
        .def("append", &PhaseList::append, "Append a new phase to the list of phase.")
        .def("data", &PhaseList::data, "Return the internal collection of Phase objects.", return_internal_ref)
        .def("empty", &PhaseList::empty, "Return true if there are no phases in the collection.")
        .def("size", &PhaseList::size, "Return the number of phases in the collection.")
        .def("species", &PhaseList::species, "Return the species that compose the phases in the collection.")
        .def("find", &PhaseList::find, "Return the index of the phase with given unique name or the number of phases if not found.")
        .def("findWithName", &PhaseList::findWithName, "Return the index of the phase with given unique name or the number of phases if not found.")
        .def("findWithSpecies", py::overload_cast<Index>(&PhaseList::findWithSpecies, py::const_), "Return the index of the phase containing the species with given name or index. If not found, return the number of phases.")
        .def("findWithSpecies", py::overload_cast<const String&>(&PhaseList::findWithSpecies, py::const_), "Return the index of the phase containing the species with given name or index. If not found, return the number of phases.")
        .def("findWithAggregateState", &PhaseList::findWithAggregateState, "Return the index of the first phase with given aggregate state of species. If not found, return the number of phases.")
        .def("findWithStateOfMatter", &PhaseList::findWithStateOfMatter, "Return the index of the first phase with given state of matter. If not found, return the number of phases.")
        .def("index", &PhaseList::index, "Return the index of the phase with given unique name. A runtime error is thrown if not found.")
        .def("indexWithName", &PhaseList::indexWithName, "Return the index of the phase with given unique name. A runtime error is thrown if not found.")
        .def("indexWithSpecies", py::overload_cast<Index>(&PhaseList::indexWithSpecies, py::const_), "Return the index of the phase containing the species with given name or index. A runtime error is thrown if not found.")
        .def("indexWithSpecies", py::overload_cast<const String&>(&PhaseList::indexWithSpecies, py::const_), "Return the index of the phase containing the species with given name or index. A runtime error is thrown if not found.")
        .def("indexWithAggregateState", &PhaseList::indexWithAggregateState, "Return the index of the first phase with given aggregate state of species. A runtime error is thrown if not found.")
        .def("indexWithStateOfMatter", &PhaseList::indexWithStateOfMatter, "Return the index of the first phase with given state of matter. A runtime error is thrown if not found.")
        .def("get", &PhaseList::get, "Return the phase with a given name.")
        .def("getWithName", &PhaseList::getWithName, "Return the phase with a given name.", return_internal_ref)
        .def("withNames", &PhaseList::withNames, "Return all phases with given names.")
        .def("withStateOfMatter", &PhaseList::withStateOfMatter, "Return all phases with given state of matter.")
        .def("withAggregateState", &PhaseList::withAggregateState, "Return all phases whose species have the given aggregate state.")
        .def("numSpeciesUntilPhase", &PhaseList::numSpeciesUntilPhase, "Return the number of species over all phases up to the one with given index.")
        .def("indicesPhasesArePure", &PhaseList::indicesPhasesArePure, "Return the indices of the phases with a single species.")
        .def("indicesPhasesAreSolution", &PhaseList::indicesPhasesAreSolution, "Return the indices of the phases with more than one species.")
        .def("indicesSpeciesInPhases", &PhaseList::indicesSpeciesInPhases, "Return the indices of the species in the given phases.")
        .def("__len__", &PhaseList::size)
        .def("__getitem__", [](const PhaseList& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const PhaseList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        .def(py::self + py::self);
        ;

    py::implicitly_convertible<Vec<Phase>, PhaseList>();
}
