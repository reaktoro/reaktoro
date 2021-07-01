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
#include <pybind11/functional.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
using namespace Reaktoro;

void exportPhase(py::module& m)
{
    py::class_<Phase>(m, "Phase")
        .def(py::init<>())
        .def("clone", &Phase::clone)
        .def("withName", &Phase::withName)
        .def("withSpecies", &Phase::withSpecies)
        .def("withStateOfMatter", &Phase::withStateOfMatter)
        .def("withActivityPropsFn", &Phase::withActivityPropsFn)
        .def("withIdealActivityPropsFn", &Phase::withIdealActivityPropsFn)
        .def("name", &Phase::name)
        .def("stateOfMatter", &Phase::stateOfMatter)
        .def("aggregateState", &Phase::aggregateState)
        .def("species", py::overload_cast<>(&Phase::species, py::const_), py::return_value_policy::reference_internal)
        .def("species", py::overload_cast<Index>(&Phase::species, py::const_), py::return_value_policy::reference_internal)
        .def("activityPropsFn", &Phase::activityPropsFn, py::return_value_policy::reference_internal)
        .def("idealActivityPropsFn", &Phase::idealActivityPropsFn, py::return_value_policy::reference_internal)
        ;
}
