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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Common/ReactionEquation.hpp>

namespace Reaktoro {

void exportReactionEquation(py::module& m)
{
    py::class_<ReactionEquation>(m, "ReactionEquation")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def(py::init<const std::map<std::string, double>&>())
        .def("numSpecies", &ReactionEquation::numSpecies)
        .def("stoichiometry", &ReactionEquation::stoichiometry)
        .def("equation", &ReactionEquation::equation, py::return_value_policy::reference_internal)
        .def("__str__", [](const ReactionEquation& x) -> std::string { return x; });
        ;
}

} // namespace Reaktoro

