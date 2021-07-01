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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Reaction.hpp>
using namespace Reaktoro;

void exportReaction(py::module& m)
{
    py::class_<Reaction>(m, "Reaction")
        .def(py::init<>())
        .def("clone", &Reaction::clone)
        .def("withName", &Reaction::withName)
        .def("withEquation", &Reaction::withEquation)
        .def("withEquilibriumConstantFn", &Reaction::withEquilibriumConstantFn)
        .def("withRateFn", &Reaction::withRateFn)
        .def("name", &Reaction::name)
        .def("equation", &Reaction::equation)
        .def("equilibriumConstantFn", &Reaction::equilibriumConstantFn)
        .def("rateFn", &Reaction::rateFn)
        ;
}
