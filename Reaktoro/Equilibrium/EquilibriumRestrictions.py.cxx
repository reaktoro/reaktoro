// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
using namespace Reaktoro;

void exportEquilibriumRestrictions(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<EquilibriumRestrictions>(m, "EquilibriumRestrictions")
        .def(py::init<const ChemicalSystem>())
        .def("cannotReact", py::overload_cast<Index>(&EquilibriumRestrictions::cannotReact))
        .def("cannotReact", py::overload_cast<String>(&EquilibriumRestrictions::cannotReact))
        .def("cannotIncrease", py::overload_cast<Index>(&EquilibriumRestrictions::cannotIncrease))
        .def("cannotIncrease", py::overload_cast<String>(&EquilibriumRestrictions::cannotIncrease))
        .def("cannotDecrease", py::overload_cast<Index>(&EquilibriumRestrictions::cannotDecrease))
        .def("cannotDecrease", py::overload_cast<String>(&EquilibriumRestrictions::cannotDecrease))
        .def("canReact", py::overload_cast<Index>(&EquilibriumRestrictions::canReact))
        .def("canReact", py::overload_cast<String>(&EquilibriumRestrictions::canReact))
        .def("canIncrease", py::overload_cast<Index>(&EquilibriumRestrictions::canIncrease))
        .def("canIncrease", py::overload_cast<String>(&EquilibriumRestrictions::canIncrease))
        .def("canDecrease", py::overload_cast<Index>(&EquilibriumRestrictions::canDecrease))
        .def("canDecrease", py::overload_cast<String>(&EquilibriumRestrictions::canDecrease))
        .def("system", &EquilibriumRestrictions::system, return_internal_ref)
        .def("speciesCannotIncrease", &EquilibriumRestrictions::speciesCannotIncrease, return_internal_ref)
        .def("speciesCannotDecrease", &EquilibriumRestrictions::speciesCannotDecrease, return_internal_ref)
        ;
}
