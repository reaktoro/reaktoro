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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Param.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data for assemblying chemical systems.
/// @ingroup Core
class Data
{
public:
    /// Construct a default Data instance with null value.
    Data();

    /// Construct a Data object as a string value.
    Data(Chars const& value);

    /// Construct a Data object as a string value.
    Data(String const& value);

    /// Construct a Data object as a number value.
    Data(double const& value);

    /// Construct a Data object as a integer value.
    Data(int const& value);

    /// Construct a Data object as a boolean value.
    Data(bool const& value);

    /// Construct a Data object as a Param object.
    Data(Param const& value);

    /// Construct a Data object as a dictionary object.
    Data(Map<String, Data> const& value);

    /// Construct a Data object as a list object.
    Data(Vec<Data> const& value);

    /// Return this data block as a string value.
    auto string() const -> String const&;

    /// Return this data block as a number value.
    auto number() const -> double const&;

    /// Return this data block as a integer value.
    auto integer() const -> int const&;

    /// Return this data block as a boolean value.
    auto boolean() const -> bool const&;

    /// Return this data block as a Param object.
    auto param() const -> Param const&;

    /// Return this data block as a dictionary object.
    auto dict() const -> Map<String, Data> const&;

    /// Return this data block as a list object.
    auto list() const -> Vec<Data> const&;

    /// Return this data block as a null value.
    auto null() const -> std::nullptr_t;

    /// Return true if this data block is a string value.
    auto isString() const -> bool;

    /// Return true if this data block is a number value.
    auto isNumber() const -> bool;

    /// Return true if this data block is a integer value.
    auto isInteger() const -> bool;

    /// Return true if this data block is a boolean value.
    auto isBoolean() const -> bool;

    /// Return true if this data block is a Param object.
    auto isParam() const -> bool;

    /// Return true if this data block is a dictionary object.
    auto isDict() const -> bool;

    /// Return true if this data block is a list object.
    auto isList() const -> bool;

    /// Return true if this data block is a null value.
    auto isNull() const -> bool;

    /// Return the child data block with given key, presuming this data block is a dictionary.
    auto operator[](String const& key) const -> Data const&;

    /// Return the child data block with given index, presuming this data block is a list.
    auto operator[](Index const& index) const -> Data const&;

    /// Add a data block to this Data object, which becomes a list if not already.
    auto add(Data const& data) -> Data&;

    /// Add a data block with given key to this Data object, which becomes a dictionary if not already.
    auto add(Chars const& key, Data const& data) -> Data&;

    /// Add a data block with given key to this Data object, which becomes a dictionary if not already.
    auto add(String const& key, Data const& data) -> Data&;

    /// Return true if a child parameter exists with given key, presuming this data block is a dictionary.
    auto exists(String const& key) const -> bool;

private:
    Any tree;
};

} // namespace Reaktoro
