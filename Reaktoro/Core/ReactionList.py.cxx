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
#include <Reaktoro/Core/ReactionList.hpp>
using namespace Reaktoro;

void exportReactionList(py::module& m)
{
    py::class_<ReactionList>(m, "ReactionList")
        .def(py::init<>())
        .def(py::init<const Vec<Reaction>&>())
        .def("append", &ReactionList::append, "Append a new reaction to the list of reactions.")
        .def("data", &ReactionList::data, "Return the internal collection of Reaction objects.", return_internal_ref)
        .def("empty", &ReactionList::empty, "Return true if there are no reactions in the collection.")
        .def("size", &ReactionList::size, "Return the number of reactions in the collection.")
        .def("find", &ReactionList::find, "Return the index of the reaction with given unique name or the number of reactions if not found.")
        .def("findWithName", &ReactionList::findWithName, "Return the index of the reaction with given unique name or the number of reactions if not found.")
        .def("index", &ReactionList::index, "Return the index of the reaction with given unique name. A runtime error is thrown if not found.")
        .def("indexWithName", &ReactionList::indexWithName, "Return the index of the reaction with given unique name. A runtime error is thrown if not found.")
        .def("get", &ReactionList::get, "Return the reaction with a given name.")
        .def("getWithName", &ReactionList::getWithName, "Return the reaction with a given name.", return_internal_ref)
        .def("withNames", &ReactionList::withNames, "Return all reactions with given names.")
        .def("__len__", &ReactionList::size)
        .def("__getitem__", [](const ReactionList& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const ReactionList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        .def(py::self + py::self);
        ;

    py::implicitly_convertible<Vec<Reaction>, ReactionList>();
}
