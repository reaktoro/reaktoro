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

#include "Table.hpp"

// C++ includes
#include <sstream>
#include <fstream>
#include <iomanip>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>

namespace Reaktoro {
namespace {

// Convenient alias for this translation unit.
using DataType = TableColumn::DataType;

/// Return true if type `T` is consistent with given column data type.
template<typename T>
auto isTypeConsistent(DataType datatype) -> bool
{
    switch(datatype)
    {
        case DataType::Float: return isFloatingPoint<T> || isInteger<T>;
        case DataType::Integer: return isInteger<T>;
        case DataType::String: return isSame<T, String>;
        case DataType::Boolean: return isSame<T, bool>;
        case DataType::Undefined: return false;
        default: return false;
    }
}

/// Return a string representation for a DataType value.
auto strColumnDataType(DataType datatype) -> String
{
    switch(datatype)
    {
        case DataType::Float: return "floating-point";
        case DataType::Integer: return "integer";
        case DataType::String: return "string";
        case DataType::Boolean: return "boolean";
        default: return "Undefined";
    }
}

/// General implementation of an append function to add a new value to a table column.
template<DataType type, typename T>
auto appendNewValue(Any& data, DataType& datatype, T const& value, Chars strvaluetype)
{
    if(!data.has_value())
    {
        data = Deque<T>();
        datatype = type;
    }

    errorifnot(isTypeConsistent<T>(type),
        "You cannot append a value of ", strvaluetype, " type to a table column that store values of ", strColumnDataType(type), " type. "
        "Make sure that after inserting the first value into a table column, the same value type is used for all subsequent inserts. "
        "Note that integer values can be stored as floating-point values in a table column of floats. No other conversion is supported.");

    auto& deque = std::any_cast<Deque<T>&>(data);
    deque.push_back(value);
}

/// Convert a TableColumn object to a vector of strings representing the column's data along its rows.
auto stringfyTableColumn(String colname, TableColumn const& column, Table::OutputOptions const& outputopts) -> Strings
{
    auto const& [delimiter, precision, scientific] = outputopts;

    Strings rows;
    rows.reserve(1 + column.rows()); // extra entry for column's name (first entry!)

    auto columnname_has_delimiter = colname.find(delimiter) != String::npos;
    rows.push_back(columnname_has_delimiter ? "\"" + colname + "\"" : colname); // e.g., delimiter is `,` and column name is `alpha,beta`; name becomes "alpha,beta" in the output

    std::ostringstream ss;
    ss << (scientific ? std::scientific : std::fixed);
    ss << std::showpoint;
    ss << std::setprecision(precision);

    for(auto i = 0; i < column.rows(); ++i)
    {
        ss.str("");
        std::visit([&](auto arg){ ss << arg; }, column[i]);
        rows.push_back(ss.str());
    }
    return rows;
}

/// Compute the maximum width of a table column already stringfied.
auto widthTableColumn(Strings const& rows) -> Index
{
    Index width = 0;
    for(auto const& row : rows)
        width = std::max(width, row.size());
    return width;
}

/// Assemble a vector with the width of each column in a table already stringfied.
auto widthsTableColumns(Vec<Strings> const& matrix) -> Vec<Index>
{
    return vectorize(matrix, RKT_LAMBDA(x, widthTableColumn(x)));
}

/// Output a Table object to a general stream object with given delimiter and numeric precision.
template<typename Stream>
auto outputTable(Stream& stream, Table const& table, Table::OutputOptions const& outputopts)
{
    // Stringfy all columns, creating a matrix of strings
    const Vec<Strings> matrix = vectorize(table.columns(), RKT_LAMBDA(pair, stringfyTableColumn(pair.first, pair.second, outputopts)));

    // Compute the maximum widths in each table column for improved output format with alignment
    const Vec<Index> widths = widthsTableColumns(matrix);

    // Define an auxiliary function to determine if a delimiter symbol is to be inserted or not
    const auto space = String(" ");
    const auto delim = [&](auto j) { return j > 0 ? outputopts.delimiter + space : space; };

    // Output the table, row by row
    for(auto i = 0; i < table.rows(); ++i) {
        for(auto j = 0; j < table.cols(); ++j)
            stream << delim(j) << std::setw(widths[j]) << (i < matrix[j].size() ? matrix[j][i] : ""); // matrix[j] is the j-th column, and matrix[j][i] is the i-th entry in the j-th column
        stream << "\n";
    }
}

} // anonymous namespace

TableColumn::TableColumn()
{}

auto TableColumn::appendFloat(double value) -> void
{
    appendNewValue<DataType::Float>(data, datatype, value, "floating-point");
    ++mrows;
}

auto TableColumn::appendInteger(long value) -> void
{
    appendNewValue<DataType::Integer>(data, datatype, value, "integer");
    ++mrows;
}

auto TableColumn::appendString(String const& value) -> void
{
    appendNewValue<DataType::String>(data, datatype, value, "string");
    ++mrows;
}

auto TableColumn::appendBoolean(bool value) -> void
{
    appendNewValue<DataType::Boolean>(data, datatype, value, "boolean");
    ++mrows;
}

auto TableColumn::dataType() const -> DataType
{
    return datatype;
}

auto TableColumn::floats() const -> Deque<double> const&
{
    errorif(datatype != DataType::Float, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to values a list of of floating-point type.");
    return std::any_cast<Deque<double> const&>(data);
}

auto TableColumn::floats() -> Deque<double>&
{
    errorif(datatype != DataType::Float, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to values a list of of floating-point type.");
    return std::any_cast<Deque<double>&>(data);
}

auto TableColumn::integers() const -> Deque<long> const&
{
    errorif(datatype != DataType::Integer, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of integer type.");
    return std::any_cast<Deque<long> const&>(data);
}

auto TableColumn::integers() -> Deque<long>&
{
    errorif(datatype != DataType::Integer, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of integer type.");
    return std::any_cast<Deque<long>&>(data);
}

auto TableColumn::strings() const -> Deque<String> const&
{
    errorif(datatype != DataType::String, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of string type.");
    return std::any_cast<Deque<String> const&>(data);
}

auto TableColumn::strings() -> Deque<String>&
{
    errorif(datatype != DataType::String, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of string type.");
    return std::any_cast<Deque<String>&>(data);
}

auto TableColumn::booleans() const -> Deque<bool> const&
{
    errorif(datatype != DataType::Boolean, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of boolean type.");
    return std::any_cast<Deque<bool> const&>(data);
}

auto TableColumn::booleans() -> Deque<bool>&
{
    errorif(datatype != DataType::Boolean, "You cannot convert a table column with values of ", strColumnDataType(datatype), " type to a list of values of boolean type.");
    return std::any_cast<Deque<bool>&>(data);
}

auto TableColumn::rows() const -> Index
{
    return mrows;
}

auto TableColumn::operator[](Index row) const -> std::variant<double, long, String, bool>
{
    errorifnot(row < mrows, "Given row index, ", row, ", is greater than number of rows in the table column, ", mrows, ".");
    switch(datatype)
    {
        case DataType::Float:     return std::any_cast<Deque<double> const&>(data)[row];
        case DataType::Integer:   return std::any_cast<Deque<long> const&>(data)[row];
        case DataType::String:    return std::any_cast<Deque<String> const&>(data)[row];
        case DataType::Boolean:   return std::any_cast<Deque<bool> const&>(data)[row];
        default: return NaN;
    }
}

Table::Table()
{

}

auto Table::columns() const -> Dict<String, TableColumn> const&
{
    return mcolumns;
}

auto Table::column(String const& columnname) const -> TableColumn const&
{
    auto it = mcolumns.find(columnname);
    errorif(it == mcolumns.end(), "There is no column named `", columnname, "` in Table object.");
    return it->second;
}

auto Table::column(String const& columnname) -> TableColumn&
{
    return mcolumns[columnname];
}

auto Table::operator[](String const& columnname) const -> Deque<double> const&
{
    return std::as_const(*this)[columnname];
}

auto Table::operator[](String const& columnname) -> Deque<double>&
{
    auto& col = column(columnname);
    auto const& datatype = col.dataType();
    errorif(datatype != DataType::Float,
        "You are using operator [] on a Table object to retrieve float values of a column that stores ", strColumnDataType(datatype), " values. "
        "Use one of the cast methods in class Table to achieve your desired result.");
    return col.floats();
}

auto Table::rows() const -> Index
{
    Index size = 0;
    for(auto const& [name, column]: mcolumns)
        size = std::max(size, column.rows());
    return size;
}

auto Table::cols() const -> Index
{
    return mcolumns.size();
}

auto Table::dump(OutputOptions const& outputopts) const -> String
{
    std::stringstream ss;
    outputTable(ss, *this, outputopts);
    return ss.str();
}

auto Table::save(String const& filepath, OutputOptions const& outputopts) const -> void
{
    std::ofstream file(filepath);
    outputTable(file, *this, outputopts);
}

Table::OutputOptions::OutputOptions()
: delimiter(""), precision(6), scientific(false)
{}

} // namespace Reaktoro
