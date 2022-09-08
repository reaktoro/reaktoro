// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/Surface.hpp>
using namespace Reaktoro;

void exportSurface(py::module& m)
{
    py::class_<Surface>(m, "Surface")
        .def(py::init<>())
        .def(py::init<String const&>())
        .def(py::init<String const&, String const&, String const&>())
        .def("clone", &Surface::clone, "Return a deep copy of this Surface object.")
        .def("withName", &Surface::withName, "Return a duplicate of this Surface object with replaced name.")
        .def("withPhases", &Surface::withPhases, "Return a duplicate of this Surface object with replaced names of the phases between which this surface exists.")
        .def("name", &Surface::name, "Return the unique name of this surface.")
        .def("phases", &Surface::phases, return_internal_ref, "Return the names of the phases between which this surface exists.")
        .def("equivalent", py::overload_cast<Surface const&>(&Surface::equivalent, py::const_), "Return true if this surface is equivalent to another with the same interface phases.")
        .def("equivalent", py::overload_cast<String const&, String const&>(&Surface::equivalent, py::const_), "Return true if this surface is equivalent to another with the same interface phases using phase names.")
        ;
}
