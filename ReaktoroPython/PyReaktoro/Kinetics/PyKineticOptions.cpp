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

#include "PyKineticSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Kinetics/KineticOptions.hpp>

namespace Reaktoro {

auto export_KineticOptions() -> void
{
    py::class_<KineticOutputOptions>("KineticOutputOptions")
        .def_readwrite("active", &KineticOutputOptions::active)
        .def_readwrite("format", &KineticOutputOptions::format)
        ;

    py::class_<KineticOptions>("KineticOptions")
        .def_readwrite("equilibrium", &KineticOptions::equilibrium)
        .def_readwrite("ode", &KineticOptions::ode)
        .def_readwrite("output", &KineticOptions::output)
        ;
}

} // namespace Reaktoro
