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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// Used to represent the data stored in a table column.
/// @see Table
class TableColumn
{
public:
    /// Used in the identification of the value type along each column in the table.
    enum class DataType
    {
        Float,     ///< The column data type indicating a double floating-point value value.
        Integer,   ///< The column data type indicating an integer value.
        String,    ///< The column data type indicating a string value.
        Boolean,   ///< The column data type indicating a boolean value.
        Undefined, ///< The column data type has not been specified yet.
    };

    /// Construct a default TableColumn object.
    TableColumn();

    /// Append a new floating-point value to the TableColumn object.
    /// @warning By appending a floating-point value, the column data type must be DataType::Float, otherwise a runtime error is thrown.
    auto appendFloat(double value) -> void;

    /// Append a new integer value to the TableColumn object.
    /// @warning By appending a integer value, the column data type must be DataType::Integer, otherwise a runtime error is thrown.
    auto appendInteger(long value) -> void;

    /// Append a new string value to the TableColumn object.
    /// @warning By appending a string value, the column data type must be DataType::String, otherwise a runtime error is thrown.
    auto appendString(String const& value) -> void;

    /// Append a new boolean value to the TableColumn object.
    /// @warning By appending a boolean value, the column data type must be DataType::Bool, otherwise a runtime error is thrown.
    auto appendBoolean(bool value) -> void;

    /// Get the data type of the column.
    auto dataType() const -> DataType;

    /// Convert this TableColumn object to a constant reference to its underlying `Deque<double>` object.
    /// @warning If the column data type is not DataType::Float, a runtime error is thrown.
    auto floats() const -> Deque<double> const&;

    /// Convert this TableColumn object to a mutable reference to its underlying `Deque<double>` object.
    /// @warning If the column data type is not DataType::Float, a runtime error is thrown.
    auto floats() -> Deque<double>&; // TODO: These methods in Table that return Deque<T>& are dangerous in which the user can change its length, and this will not be reflected in mrows. Replace this by std::span when migrating to C++20.

    /// Convert this TableColumn object to a constant reference to its underlying `Deque<long>` object.
    /// @warning If the column data type is not DataType::Integer, a runtime error is thrown.
    auto integers() const -> Deque<long> const&;

    /// Convert this TableColumn object to a mutable reference to its underlying `Deque<long>` object.
    /// @warning If the column data type is not DataType::Integer, a runtime error is thrown.
    auto integers() -> Deque<long>&;

    /// Convert this TableColumn object to a constant reference to its underlying `Deque<String>` object.
    /// @warning If the column data type is not DataType::String, a runtime error is thrown.
    auto strings() const -> Deque<String> const&;

    /// Convert this TableColumn object to a mutable reference to its underlying `Deque<String>` object.
    /// @warning If the column data type is not DataType::String, a runtime error is thrown.
    auto strings() -> Deque<String>&;

    /// Convert this TableColumn object to a constant reference to its underlying `Deque<bool>` object.
    /// @warning If the column data type is not DataType::Bool, a runtime error is thrown.
    auto booleans() const -> Deque<bool> const&;

    /// Convert this TableColumn object to a mutable reference to its underlying `Deque<bool>` object.
    /// @warning If the column data type is not DataType::Bool, a runtime error is thrown.
    auto booleans() -> Deque<bool>&;

    /// Get the number of rows in the column.
    auto rows() const -> Index;

    /// Get the value stored in the column at the given row index.
    auto operator[](Index row) const -> std::variant<double, long, String, bool>;

    /// Append a new value to the TableColumn object.
    /// This method is tolerant of inserting an integer value into a table column of floats. In this
    /// case, the integer value is converted to a float.
    template<typename T>
    auto append(T const& value) -> void
    {
        if constexpr(isSame<T, bool>)
            appendBoolean(value);
        else if constexpr(isFloatingPoint<T>)
            appendFloat(value);
        else if constexpr(isInteger<T>) // keep integer check after booleans, as isInteger<bool> is true!
            if(datatype == DataType::Float)
                appendFloat(value); // allow integers to be cast to float and stored in a float column!
            else appendInteger(value);
        else if constexpr(isSame<T, String>)
            appendString(value);
        else errorif(true, "You cannot append this value with an unsupported type to a table column.");
    }

    /// Append a new string value to the TableColumn object.
    auto append(Chars value) -> void
    {
        appendString(value);
    }

    /// Append a new value to the TableColumn object.
    /// This operator exists for convenience and provides an operation almost identical to @ref
    /// append. The difference is that this operator treats integers as floats. This is to avoid
    /// accidental first insertion of an integer value, but subsequent inserted values are floats
    /// (especially from the Python side; see example below). If you want a table column of
    /// integers, use the @ref appendInteger or @ref append method with an integer value.
    template<typename T>
    auto operator<<(T const& value) -> TableColumn&
    {
        if constexpr(isSame<T, bool>)
            appendBoolean(value);
        else if constexpr(isFloatingPoint<T> || isInteger<T>) // keep integer check after booleans, as isInteger<bool> is true!
            appendFloat(value); // note that integers are cast to float with this operator!
        else if constexpr(isSame<T, String>)
            appendString(value);
        else errorif(true, "You cannot append this value with an unsupported type to a table column.");
        return *this;
    }

    /// Append a new string value to the TableColumn object.
    auto operator<<(Chars value) -> TableColumn&
    {
        appendString(value);
        return *this;
    }

    /// Cast this TableColumn object to a mutable reference to a list of values with type compatible with given one.
    template<typename T>
    auto cast() -> Deque<T>&
    {
        if constexpr(isSame<T, bool>)
            return booleans();
        else if constexpr(isFloatingPoint<T>)
            return floats();
        else if constexpr(isInteger<T>) // keep integer check after booleans, as isInteger<bool> is true!
            return integers();
        else if constexpr(isSame<T, String> || isSame<T, Chars>)
            return strings();
        else errorif(true, "You cannot cast this table column to a list of values with an unsupported type.");
    }

    /// Cast this TableColumn object to a constant reference to a list of values with type compatible with given one.
    template<typename T>
    auto cast() const -> Deque<T> const&
    {
        return std::as_const(*this).cast<T>();
    }

private:
    /// The values stored in this table column (e.g., `Deque<double>`, `Deque<long>`, `Deque<String>`, `Deque<bool>`).
    Any data;

    /// The number of rows in the column.
    Index mrows = 0;

    /// The type of values stored in this table column.
    DataType datatype = DataType::Undefined;
};

/// Used to store computed data in columns.
class Table
{
public:
    /// Construct a default Table object.
    Table();

    /// Get the columns in the Table object.
    auto columns() const -> Dict<String, TableColumn> const&;

    /// Get a constant reference to a column in the table with given name.
    auto column(String const& columnname) const -> TableColumn const&;

    /// Get a mutable reference to a column in the table with given name.
    auto column(String const& columnname) -> TableColumn&;

    /// Get a constant reference to a column in the table with given name.
    auto operator[](String const& columnname) const -> Deque<double> const&;

    /// Get a mutable reference to a column in the table with given name.
    auto operator[](String const& columnname) -> Deque<double>&;

    /// Get the number of rows in the table (i.e., the length of the longest column in the table).
    auto rows() const -> Index;

    /// Get the number of columns in the table.
    auto cols() const -> Index;

    /// Used to specify formatting options when printing or outputting the Table object.
    struct OutputOptions
    {
        /// Construct a default OutputOptions object.
        OutputOptions();

        /// The symbol used to separate column values on a table row (defaults to `""`).
        String delimiter;

        /// The precision used when printing floating-point values (defaults to 6).
        int precision;

        /// The boolean flag indicating if floating-point values should be printed in scientific notation (defaults to `false`).
        bool scientific;
    };

    /// Assemble a string representation of the Table object.
    /// @param outputopts The formatting options for the output operation.
    auto dump(OutputOptions const& outputopts = {}) const -> String;

    /// Save the Table object to a file.
    /// @param filepath The path to the file that will be created, including its file name (e.g., `table.txt`, `/home/username/downloads/table.txt`).
    /// @param outputopts The formatting options for the output operation.
    /// @warning Ensure that the path given exists; no directories are created in this method call.
    auto save(String const& filepath, OutputOptions const& outputopts = {}) const -> void;

private:
    /// The named columns and their stored values in the table.
    Dict<String, TableColumn> mcolumns;
};

} // namespace Reaktoro
