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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/PhaseList.hpp>
using namespace Reaktoro;

template<typename PhaseListType>
auto createCommonPhaseList(py::module& m, const char* name)
{
    return
    py::class_<PhaseListType>(m, name)
        .def("data", &PhaseListType::data)
        .def("size", &PhaseListType::size)
        .def("find", &PhaseListType::find)
        .def("findWithName", &PhaseListType::findWithName)
        .def("findWithSpecies", py::overload_cast<Index>(&PhaseListType::findWithSpecies, py::const_))
        .def("findWithSpecies", py::overload_cast<const String&>(&PhaseListType::findWithSpecies, py::const_))
        .def("findWithAggregateState", &PhaseListType::findWithAggregateState)
        .def("findWithStateOfMatter", &PhaseListType::findWithStateOfMatter)
        .def("index", &PhaseListType::index)
        .def("indexWithName", &PhaseListType::indexWithName)
        .def("indexWithSpecies", py::overload_cast<Index>(&PhaseListType::indexWithSpecies, py::const_))
        .def("indexWithSpecies", py::overload_cast<const String&>(&PhaseListType::indexWithSpecies, py::const_))
        .def("indexWithAggregateState", &PhaseListType::indexWithAggregateState)
        .def("indexWithStateOfMatter", &PhaseListType::indexWithStateOfMatter)
        .def("withNames", &PhaseListType::withNames)
        .def("withStateOfMatter", &PhaseListType::withStateOfMatter)
        .def("withAggregateState", &PhaseListType::withAggregateState)
        .def("numSpeciesUntilPhase", &PhaseListType::numSpeciesUntilPhase)
        .def("__getitem__", [](const PhaseListType& self, Index i) { return self[i]; }, py::return_value_policy::reference_internal)
        .def("__iter__", [](const PhaseListType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()); // keep object alive while iterator exists;
        ;
}

void exportPhaseList(py::module& m)
{
    createCommonPhaseList<PhaseList>(m, "PhaseList")
        .def(py::init<>())
        .def(py::init<const Vec<Phase>&>())
        .def("append", &PhaseList::append)
        ;

    createCommonPhaseList<PhaseListConstRef>(m, "PhaseListConstRef");
}
