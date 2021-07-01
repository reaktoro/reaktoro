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
        .def("cannotIncreaseAbove", py::overload_cast<Index, double, String>(&EquilibriumRestrictions::cannotIncreaseAbove))
        .def("cannotIncreaseAbove", py::overload_cast<String, double, String>(&EquilibriumRestrictions::cannotIncreaseAbove))
        .def("cannotDecrease", py::overload_cast<Index>(&EquilibriumRestrictions::cannotDecrease))
        .def("cannotDecrease", py::overload_cast<String>(&EquilibriumRestrictions::cannotDecrease))
        .def("cannotDecreaseBelow", py::overload_cast<Index, double, String>(&EquilibriumRestrictions::cannotDecreaseBelow))
        .def("cannotDecreaseBelow", py::overload_cast<String, double, String>(&EquilibriumRestrictions::cannotDecreaseBelow))
        .def("canReactFreely", py::overload_cast<Index>(&EquilibriumRestrictions::canReactFreely))
        .def("canReactFreely", py::overload_cast<String>(&EquilibriumRestrictions::canReactFreely))
        .def("canIncreaseFreely", py::overload_cast<Index>(&EquilibriumRestrictions::canIncreaseFreely))
        .def("canIncreaseFreely", py::overload_cast<String>(&EquilibriumRestrictions::canIncreaseFreely))
        .def("canDecreaseFreely", py::overload_cast<Index>(&EquilibriumRestrictions::canDecreaseFreely))
        .def("canDecreaseFreely", py::overload_cast<String>(&EquilibriumRestrictions::canDecreaseFreely))
        .def("system", &EquilibriumRestrictions::system, return_internal_ref)
        .def("speciesCannotIncrease", &EquilibriumRestrictions::speciesCannotIncrease, return_internal_ref)
        .def("speciesCannotDecrease", &EquilibriumRestrictions::speciesCannotDecrease, return_internal_ref)
        .def("speciesCannotIncreaseAbove", &EquilibriumRestrictions::speciesCannotIncreaseAbove, return_internal_ref)
        .def("speciesCannotDecreaseBelow", &EquilibriumRestrictions::speciesCannotDecreaseBelow, return_internal_ref)
        ;
}
