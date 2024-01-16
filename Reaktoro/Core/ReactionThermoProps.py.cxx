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
#include <Reaktoro/Core/ReactionThermoProps.hpp>
using namespace Reaktoro;

void exportReactionThermoProps(py::module& m)
{
    py::class_<ReactionThermoProps>(m, "ReactionThermoProps")
        .def(py::init<>())
        .def_readwrite("T"   , &ReactionThermoProps::T)
        .def_readwrite("P"   , &ReactionThermoProps::P)
        .def_readwrite("lgK" , &ReactionThermoProps::lgK)
        .def_readwrite("dG0" , &ReactionThermoProps::dG0)
        .def_readwrite("dH0" , &ReactionThermoProps::dH0)
        .def_readwrite("dV0" , &ReactionThermoProps::dV0)
        .def_readwrite("dVT0", &ReactionThermoProps::dVT0)
        .def_readwrite("dVP0", &ReactionThermoProps::dVP0)
        .def_readwrite("dCp0", &ReactionThermoProps::dCp0)
        .def_readwrite("dCv0", &ReactionThermoProps::dCv0)
        .def_readwrite("dU0" , &ReactionThermoProps::dU0)
        .def_readwrite("dS0" , &ReactionThermoProps::dS0)
        .def_readwrite("dA0" , &ReactionThermoProps::dA0)
        .def("__repr__", [](const ReactionThermoProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
