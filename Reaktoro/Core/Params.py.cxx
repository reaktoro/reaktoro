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
#include <Reaktoro/Core/Params.hpp>
using namespace Reaktoro;

void exportParams(py::module& m)
{
    py::class_<Params>(m, "Params")
        .def(py::init<>(), "Construct a default Params object.")
        .def(py::init<Data const&>(), "Construct a Params object with given parameters as Data object.")
        .def_static("fromFile", &Params::fromFile, "Return a Params object constructed with a given local file.")
        .def_static("fromEmbeddedFile", &Params::fromEmbeddedFile, "Return a Params object constructed with a given embedded file.")
        .def_static("fromStringYAML", &Params::fromStringYAML, "Return a Params object constructed with given parameters text contents in YAML format.")
        .def_static("fromStringJSON", &Params::fromStringJSON, "Return a Params object constructed with given parameters text contents in JSON format.")
        .def_static("embedded", &Params::embedded, "Return parameters with the given file path within the virtual directory of embedded resources.")
        .def_static("local", &Params::local, "Return parameters with the given local file path.")
        .def("data", &Params::data, return_internal_ref, "Return the underlying Data object of this Params object.")
        .def("append", py::overload_cast<Params const&>(&Params::append), return_internal_ref, "Append another Params object into this.")
        .def("append", py::overload_cast<Data const&>(&Params::append), return_internal_ref, "Append a Data object containing parameters into this Params object.")
        .def("__iadd__", [](Params& self, Params const& other) { return self.append(other); }, return_internal_ref, "Append another Params object into this.")
        .def("__getitem__", [](Params& self, String const& name) { return self[name]; }, return_internal_ref, "Return the set of parameters with given name.")
        ;
}
