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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ReactionProps.hpp>
using namespace Reaktoro;

void exportReactionProps(py::module& m)
{
    py::class_<ReactionProps>(m, "ReactionProps")
        .def(py::init<>())
        .def_readwrite("T"   , &ReactionProps::T)
        .def_readwrite("P"   , &ReactionProps::P)
        .def_readwrite("lgK" , &ReactionProps::lgK)
        .def_readwrite("dG0" , &ReactionProps::dG0)
        .def_readwrite("dH0" , &ReactionProps::dH0)
        .def_readwrite("dV0" , &ReactionProps::dV0)
        .def_readwrite("dVT0", &ReactionProps::dVT0)
        .def_readwrite("dVP0", &ReactionProps::dVP0)
        .def_readwrite("dCp0", &ReactionProps::dCp0)
        .def_readwrite("dCv0", &ReactionProps::dCv0)
        .def_readwrite("dU0" , &ReactionProps::dU0)
        .def_readwrite("dS0" , &ReactionProps::dS0)
        .def_readwrite("dA0" , &ReactionProps::dA0)
        .def("__repr__", [](const ReactionProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
