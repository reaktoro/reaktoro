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
        .def("withActivityModel", &Phase::withActivityModel)
        .def("withIdealActivityModel", &Phase::withIdealActivityModel)
        .def("name", &Phase::name)
        .def("stateOfMatter", &Phase::stateOfMatter)
        .def("aggregateState", &Phase::aggregateState)
        .def("elements", &Phase::elements, py::return_value_policy::reference_internal)
        .def("element", &Phase::element, py::return_value_policy::reference_internal)
        .def("species", py::overload_cast<>(&Phase::species, py::const_), py::return_value_policy::reference_internal)
        .def("species", py::overload_cast<Index>(&Phase::species, py::const_), py::return_value_policy::reference_internal)
        .def("activityModel", &Phase::activityModel, py::return_value_policy::reference_internal)
        .def("idealActivityModel", &Phase::idealActivityModel, py::return_value_policy::reference_internal)
        ;
}
