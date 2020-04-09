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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Kinetics/KineticOptions.hpp>

namespace Reaktoro {

void exportKineticOptions(py::module& m)
{
    py::class_<KineticOutputOptions>(m, "KineticOutputOptions")
        .def_readwrite("active", &KineticOutputOptions::active)
        .def_readwrite("format", &KineticOutputOptions::format)
        ;

    py::class_<KineticOptions>(m, "KineticOptions")
        .def(py::init<>())
        .def_readwrite("equilibrium", &KineticOptions::equilibrium)
        .def_readwrite("ode", &KineticOptions::ode)
        .def_readwrite("output", &KineticOptions::output)
        ;
}

} // namespace Reaktoro
