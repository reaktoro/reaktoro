// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/YAML.hpp>
using namespace Reaktoro;

void exportYAML(py::module& m)
{
    py::class_<yaml>(m, "yaml")
        .def(py::init<>())
        .def_static("parse", py::overload_cast<const char*>(&yaml::parse))
        .def_static("parse", py::overload_cast<const std::string&>(&yaml::parse))
        .def_static("parse", py::overload_cast<std::istream&>(&yaml::parse))
        .def("__getitem__", [](const yaml& self, const std::string& key) { return self[key]; } )
        .def("__getitem__", [](yaml& self, const std::string& key) { return self[key]; } )
        .def("__repr__", &yaml::repr)
        ;
}
