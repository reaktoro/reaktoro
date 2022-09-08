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
#include <Reaktoro/Core/Data.hpp>
using namespace Reaktoro;

void exportData(py::module& m)
{
    py::class_<Data>(m, "Data")
        .def(py::init<>(), "Construct a default Data instance with null value.")
        .def(py::init<Chars const&>(), "Construct a Data object as a string value.")
        .def(py::init<String const&>(), "Construct a Data object as a string value.")
        .def(py::init<double const&>(), "Construct a Data object as a number value.")
        .def(py::init<int const&>(), "Construct a Data object as a integer value.")
        .def(py::init<bool const&>(), "Construct a Data object as a boolean value.")
        .def(py::init<Param const&>(), "Construct a Data object as a Param object.")
        .def(py::init<Map<String, Data> const&>(), "Construct a Data object as a dictionary object.")
        .def(py::init<Vec<Data> const&>(), "Construct a Data object as a list object.")
        .def("string", &Data::string, return_internal_ref, "Return this data block as a string value.")
        .def("number", &Data::number, return_internal_ref, "Return this data block as a number value.")
        .def("integer", &Data::integer, return_internal_ref, "Return this data block as a integer value.")
        .def("boolean", [](Data const& self) { return self.integer() != 0; }, "Return this data block as a boolean value.")
        .def("param", &Data::param, return_internal_ref, "Return this data block as a Param object.")
        .def("dict", &Data::dict, return_internal_ref, "Return this data block as a dictionary object.")
        .def("list", &Data::list, return_internal_ref, "Return this data block as a list object.")
        .def("null", &Data::null, "Return this data block as a null value.")
        .def("isString", &Data::isString, "Return true if this data block is a string value.")
        .def("isNumber", &Data::isNumber, "Return true if this data block is a number value.")
        .def("isInteger", &Data::isInteger, "Return true if this data block is a integer value.")
        .def("isBoolean", [](Data const& self) { return self.isInteger(); }, "Return true if this data block is a boolean value.")
        .def("isParam", &Data::isParam, "Return true if this data block is a Param object.")
        .def("isDict", &Data::isDict, "Return true if this data block is a dictionary object.")
        .def("isList", &Data::isList, "Return true if this data block is a list object.")
        .def("isNull", &Data::isNull, "Return true if this data block is a null value.")
        .def("__getitem__", [](Data const& self, String const& key) { return self[key]; }, "Return the child data block with given key, presuming this data block is a dictionary.")
        .def("__getitem__", [](Data const& self, Index const& index) { return self[index]; }, "Return the child data block with given index, presuming this data block is a list.")
        .def("with", &Data::with, "Return the child data block whose `attribute` has a given `value`, presuming this data block is a list.")
        .def("add", py::overload_cast<Data const&>(&Data::add), "Add a data block to this Data object, which becomes a list if not already.")
        .def("add", py::overload_cast<Chars const&, Data const&>(&Data::add), "Add a data block with given key to this Data object, which becomes a dictionary if not already.")
        .def("add", py::overload_cast<String const&, Data const&>(&Data::add), "Add a data block with given key to this Data object, which becomes a dictionary if not already.")
        .def("exists", &Data::exists, "Return true if a child parameter exists with given key, presuming this data block is a dictionary.")
        ;

    py::implicitly_convertible<Chars, Data>();
    py::implicitly_convertible<String, Data>();
    py::implicitly_convertible<double, Data>();
    py::implicitly_convertible<int, Data>();
    py::implicitly_convertible<bool, Data>();
    py::implicitly_convertible<Param, Data>();
    py::implicitly_convertible<Map<String, Data>, Data>();
    py::implicitly_convertible<Vec<Data>, Data>();
}
