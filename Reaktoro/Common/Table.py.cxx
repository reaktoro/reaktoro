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
#include <Reaktoro/Common/Table.hpp>
using namespace Reaktoro;

template<typename T>
auto __lshift__(TableColumn& self, T const& value) -> TableColumn&
{
    return self << value;
}

void exportTable(py::module& m)
{
    // Create a module object for TableColumn
    auto mTableColumn = py::class_<TableColumn>(m, "TableColumn");

    // Register the TableColumn.DataType nested type first to avoid potential "type not registered yet" errors.
    py::enum_<TableColumn::DataType>(mTableColumn, "DataType")
        .value("Float", TableColumn::DataType::Float, "The column data type indicating a double floating-point value value.")
        .value("Integer", TableColumn::DataType::Integer, "The column data type indicating an integer value.")
        .value("String", TableColumn::DataType::String, "The column data type indicating a string value.")
        .value("Boolean", TableColumn::DataType::Boolean, "The column data type indicating a boolean value.")
        .value("Undefined", TableColumn::DataType::Undefined, "The column data type has not been specified yet.")
        ;

    // Finalize registration of remaining TableColumn methods and attributes.
    mTableColumn
        .def(py::init<>())
        .def("appendFloat", &TableColumn::appendFloat, "Append a new floating-point value to the TableColumn object.")
        .def("appendInteger", &TableColumn::appendInteger, "Append a new integer value to the TableColumn object.")
        .def("appendString", &TableColumn::appendString, "Append a new string value to the TableColumn object.")
        .def("appendBoolean", &TableColumn::appendBoolean, "Append a new boolean value to the TableColumn object.")
        .def("dataType", &TableColumn::dataType, "Get the data type of the column.")
        .def("floats", py::overload_cast<>(&TableColumn::floats, py::const_), return_internal_ref, "Convert this TableColumn object to a constant reference to its underlying `Deque<double>` object.")
        .def("floats", py::overload_cast<>(&TableColumn::floats), return_internal_ref, "Convert this TableColumn object to a mutable reference to its underlying `Deque<double>` object.")
        .def("integers", py::overload_cast<>(&TableColumn::integers, py::const_), return_internal_ref, "Convert this TableColumn object to a constant reference to its underlying `Deque<long>` object.")
        .def("integers", py::overload_cast<>(&TableColumn::integers), return_internal_ref, "Convert this TableColumn object to a mutable reference to its underlying `Deque<long>` object.")
        .def("strings", py::overload_cast<>(&TableColumn::strings, py::const_), return_internal_ref, "Convert this TableColumn object to a constant reference to its underlying `Deque<String>` object.")
        .def("strings", py::overload_cast<>(&TableColumn::strings), return_internal_ref, "Convert this TableColumn object to a mutable reference to its underlying `Deque<String>` object.")
        .def("booleans", py::overload_cast<>(&TableColumn::booleans, py::const_), return_internal_ref, "Convert this TableColumn object to a constant reference to its underlying `Deque<bool>` object.")
        .def("booleans", py::overload_cast<>(&TableColumn::booleans), return_internal_ref, "Convert this TableColumn object to a mutable reference to its underlying `Deque<bool>` object.")
        .def("rows", &TableColumn::rows, "Get the number of rows in the column.")
        .def("__getitem__", [](TableColumn const& self, int irow) { return self[irow]; } )
        .def("append", [](TableColumn& self, bool value) { self.append(value); }, "Append a new value to the TableColumn object.")
        .def("append", [](TableColumn& self, double value) { self.append(value); }, "Append a new value to the TableColumn object.")
        .def("append", [](TableColumn& self, long value) { self.append(value); }, "Append a new value to the TableColumn object.")
        .def("append", [](TableColumn& self, String const& value) { self.append(value); }, "Append a new value to the TableColumn object.")
        .def("__lshift__", __lshift__<bool>, return_internal_ref, "Append a new value to the TableColumn object.")
        .def("__lshift__", __lshift__<double>, return_internal_ref, "Append a new value to the TableColumn object.")
        .def("__lshift__", __lshift__<long>, return_internal_ref, "Append a new value to the TableColumn object.")
        .def("__lshift__", __lshift__<String>, return_internal_ref, "Append a new value to the TableColumn object.")
        ;

    // Create a module object for TableColumn
    auto mTable = py::class_<Table>(m, "Table");

    // Register the Table.OutputOptions nested type first to avoid potential "type not registered yet" errors.
    py::class_<Table::OutputOptions>(mTable, "OutputOptions")
        .def(py::init<>())
        .def_readwrite("delimiter", &Table::OutputOptions::delimiter, "The symbol used to separate column values on a table row (defaults to `;`).")
        .def_readwrite("precision", &Table::OutputOptions::precision, "The precision used when printing floating-point values (defaults to 6).")
        .def_readwrite("scientific", &Table::OutputOptions::scientific, "The boolean flag indicating if floating-point values should be printed in scientific notation (defaults to `false`).")
        .def_readwrite("fixed", &Table::OutputOptions::fixed, "The boolean flag indicating if floating-point values should be printed in fixed notation (defaults to `false`).")
        ;

    // Finalize registration of remaining Table methods and attributes.
    mTable
        .def(py::init<>())
        .def("columns", &Table::columns, "Get the columns in the Table object.")
        .def("column", py::overload_cast<String const&>(&Table::column), return_internal_ref, "Get a mutable reference to a column in the table with given name.")
        .def("__getitem__", [](Table& self, String const& colname) { return self[colname]; }, return_internal_ref, "Get a mutable reference to a column in the table with given name.")
        .def("rows", &Table::rows, "Get the number of rows in the table (i.e., the length of the longest column in the table).")
        .def("cols", &Table::cols, "Get the number of columns in the table.")
        .def("dump", &Table::dump, "Assemble a string representation of the Table object.", "outputopts"_a = Table::OutputOptions())
        .def("save", &Table::save, "Save the Table object to a file.", "filepath"_a, "outputopts"_a = Table::OutputOptions())
        .def("__str__", [](Table const& self) { return self.dump(); })
        ;
}
