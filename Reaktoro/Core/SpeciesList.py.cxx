// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/SpeciesList.hpp>
using namespace Reaktoro;

template<typename SpeciesListType>
auto createCommonSpeciesList(py::module& m, const char* name)
{
    return
    py::class_<SpeciesListType>(m, name)
        .def("data", &SpeciesListType::data)
        .def("size", &SpeciesListType::size)
        .def("find", &SpeciesListType::find)
        .def("findWithName", &SpeciesListType::findWithName)
        .def("findWithFormula", &SpeciesListType::findWithFormula)
        .def("findWithSubstance", &SpeciesListType::findWithSubstance)
        .def("index", &SpeciesListType::index)
        .def("indexWithName", &SpeciesListType::indexWithName)
        .def("indexWithFormula", &SpeciesListType::indexWithFormula)
        .def("indexWithSubstance", &SpeciesListType::indexWithSubstance)
        .def("get", &SpeciesListType::get, py::return_value_policy::reference_internal)
        .def("getWithName", &SpeciesListType::getWithName, py::return_value_policy::reference_internal)
        .def("getWithFormula", &SpeciesListType::getWithFormula, py::return_value_policy::reference_internal)
        .def("getWithSubstance", &SpeciesListType::getWithSubstance, py::return_value_policy::reference_internal)
        .def("withNames", &SpeciesListType::withNames)
        .def("withFormulas", &SpeciesListType::withFormulas)
        .def("withSubstances", &SpeciesListType::withSubstances)
        .def("withAggregateState", &SpeciesListType::withAggregateState)
        .def("withTag", &SpeciesListType::withTag)
        .def("withoutTag", &SpeciesListType::withoutTag)
        .def("withTags", &SpeciesListType::withTags)
        .def("withoutTags", &SpeciesListType::withoutTags)
        .def("withElements", &SpeciesListType::withElements)
        .def("withElementsOf", &SpeciesListType::withElementsOf)
        .def("__getitem__", [](const SpeciesListType& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const SpeciesListType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;
}

void exportSpeciesList(py::module& m)
{
    createCommonSpeciesList<SpeciesList>(m, "SpeciesList")
        .def(py::init<>())
        .def(py::init<const StringList&>())
        .def(py::init<const Vec<Species>&>())
        .def("append", &SpeciesList::append)
        ;

    createCommonSpeciesList<SpeciesListConstRef>(m, "SpeciesListConstRef");

    py::implicitly_convertible<Vec<Species>, SpeciesList>();
    py::implicitly_convertible<Vec<Species>, SpeciesListConstRef>();
    py::implicitly_convertible<SpeciesList, SpeciesListConstRef>();
}
