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

// Third-party includes
#include <nlohmann/json_fwd.hpp>

// Third-party forward declarations
namespace YAML { class Node; }

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Core/Param.hpp>

namespace Reaktoro {

// ================================================================================================
// TODO: Implement Data::with(key_name, key_value) to find entries in lists. Example:
//
// params["Species"].with("Name", "H2O(aq)")["StandardThermoModel"]["HKF"]["Gf"]
//
// Implement also Data::get(string_path_to_parameter) such as:
//
// params.get("Species:Name=H2O(aq):StandardThermoModel:HKF:Gf")
// params.get("Species:Formula=Ca++:StandardThermoModel:HKF:Gf")
// ================================================================================================

/// The class used to store and retrieve data for assemblying chemical systems.
/// @ingroup Core
class Data
{
public:
    /// Construct a default Data instance with null value.
    Data();

    /// Construct a Data object with given boolean value.
    Data(bool const& value);

    /// Construct a Data object with given string value.
    Data(Chars const& value);

    /// Construct a Data object with given string value.
    Data(String const& value);

    /// Construct a Data object with given Param object.
    Data(Param const& value);

    /// Construct a Data object with given numeric value (integers and floating-point values), which becomes a Param object.
    template<typename Number, EnableIf<isNumeric<Number>>...>
    Data(Number const& value) : Data(Param(value)) {}

    /// Construct a Data object with given dictionary object.
    Data(Map<String, Data> const& value);

    /// Construct a Data object with given list object.
    Data(Vec<Data> const& value);

    /// Construct a Data object with given yaml object.
    Data(YAML::Node const& obj);

    /// Construct a Data object with given json object.
    Data(nlohmann::json const& obj);

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(Chars const& input) -> Data;

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(String const& input) -> Data;

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(std::istream& input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(Chars const& input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(String const& input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(std::istream& input) -> Data;

    /// Return this data block as a string value.
    auto string() const -> String const&;

    /// Return this data block as a number value.
    auto number() const -> double;

    /// Return this data block as a integer value.
    auto integer() const -> int;

    /// Return this data block as a boolean value.
    auto boolean() const -> bool;

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
