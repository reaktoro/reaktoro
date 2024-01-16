// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Reaction.hpp>
using namespace Reaktoro;

void exportReaction(py::module& m)
{
    py::class_<Reaction>(m, "Reaction")
        .def(py::init<>())
        .def("clone", &Reaction::clone)
        .def("withName", &Reaction::withName)
        .def("withEquation", &Reaction::withEquation)
        .def("withRateModel", &Reaction::withRateModel)
        .def("name", &Reaction::name)
        .def("equation", &Reaction::equation)
        .def("rateModel", &Reaction::rateModel)
        .def("props", py::overload_cast<real, real>(&Reaction::props, py::const_))
        .def("props", py::overload_cast<real, Chars, real, Chars>(&Reaction::props, py::const_))
        .def("rate", &Reaction::rate)
        ;
}
