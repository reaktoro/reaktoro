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
        .def("get", &SpeciesList::get, py::return_value_policy::reference_internal)
        .def("getWithName", &SpeciesList::getWithName, py::return_value_policy::reference_internal)
        .def("getWithFormula", &SpeciesList::getWithFormula, py::return_value_policy::reference_internal)
        .def("getWithSubstance", &SpeciesList::getWithSubstance, py::return_value_policy::reference_internal)
        .def("withNames", &SpeciesList::withNames)
        .def("withFormulas", &SpeciesList::withFormulas)
        .def("withSubstances", &SpeciesList::withSubstances)
        .def("withAggregateState", &SpeciesList::withAggregateState)
        .def("withTag", &SpeciesList::withTag)
        .def("withoutTag", &SpeciesList::withoutTag)
        .def("withTags", &SpeciesList::withTags)
        .def("withoutTags", &SpeciesList::withoutTags)
        .def("withElements", &SpeciesList::withElements)
        .def("withElementsOf", &SpeciesList::withElementsOf)
        .def("__getitem__", [](const SpeciesList& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const SpeciesList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;

    py::implicitly_convertible<Vec<Species>, SpeciesList>();
}