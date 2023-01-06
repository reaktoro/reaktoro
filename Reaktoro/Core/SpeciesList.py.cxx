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
#include <Reaktoro/Core/SpeciesList.hpp>
using namespace Reaktoro;

void exportSpeciesList(py::module& m)
{
    py::class_<SpeciesList>(m, "SpeciesList")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Vec<Species>&>())
        .def("append", &SpeciesList::append)
        .def("data", &SpeciesList::data)
        .def("empty", &SpeciesList::empty)
        .def("size", &SpeciesList::size)
        .def("elements", &SpeciesList::elements)
        .def("find", &SpeciesList::find)
        .def("findWithName", &SpeciesList::findWithName)
        .def("findWithFormula", &SpeciesList::findWithFormula)
        .def("findWithSubstance", &SpeciesList::findWithSubstance)
        .def("index", &SpeciesList::index)
        .def("indexWithName", &SpeciesList::indexWithName)
        .def("indexWithFormula", &SpeciesList::indexWithFormula)
        .def("indexWithSubstance", &SpeciesList::indexWithSubstance)
        .def("front", &SpeciesList::front, return_internal_ref)
        .def("back", &SpeciesList::back, return_internal_ref)
        .def("get", &SpeciesList::get, return_internal_ref)
        .def("getWithName", &SpeciesList::getWithName, return_internal_ref)
        .def("getWithFormula", &SpeciesList::getWithFormula, return_internal_ref)
        .def("getWithSubstance", &SpeciesList::getWithSubstance, return_internal_ref)
        .def("withNames", &SpeciesList::withNames)
        .def("withFormulas", &SpeciesList::withFormulas)
        .def("withSubstances", &SpeciesList::withSubstances)
        .def("withAggregateState", &SpeciesList::withAggregateState)
        .def("withCharge", &SpeciesList::withCharge)
        .def("withTag", &SpeciesList::withTag)
        .def("withoutTag", &SpeciesList::withoutTag)
        .def("withTags", &SpeciesList::withTags)
        .def("withoutTags", &SpeciesList::withoutTags)
        .def("withElements", &SpeciesList::withElements)
        .def("withElementsOf", &SpeciesList::withElementsOf)
        .def("__len__", &SpeciesList::size)
        .def("__getitem__", [](const SpeciesList& self, Index i) { return self[i]; }, return_internal_ref)
        .def("__iter__", [](const SpeciesList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        .def(py::self + py::self);
        ;

    py::implicitly_convertible<Vec<Species>, SpeciesList>();
}
