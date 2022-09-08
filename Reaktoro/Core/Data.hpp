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

    /// Construct a Data object with given yaml object.
    Data(YAML::Node const& obj);

    /// Construct a Data object with given json object.
    Data(nlohmann::json const& obj);

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(Chars input) -> Data;

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(String const& input) -> Data;

    /// Return a Data object by parsing an YAML document.
    static auto fromYaml(std::istream& input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(Chars input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(String const& input) -> Data;

    /// Return a Data object by parsing a JSON document.
    static auto fromJson(std::istream& input) -> Data;

    /// Return this Data object as a string.
    auto asString() const -> String const&;

    /// Return this Data object as a float number.
    auto asFloat() const -> double;

    /// Return this Data object as an integer number.
    auto asInteger() const -> int;

    /// Return this Data object as a boolean value.
    auto asBoolean() const -> bool;

    /// Return this Data object as a real number.
    auto asReal() const -> real const&;

    /// Return this Data object as a Param object.
    auto asParam() const -> Param const&;

    /// Return this Data object as a dictionary object.
    auto asDict() const -> Map<String, Data> const&;

    /// Return this Data object as a list object.
    auto asList() const -> Vec<Data> const&;

    /// Return true if this Data object is a string.
    auto isString() const -> bool;

    /// Return true if this Data object is a float number.
    auto isFloat() const -> bool;

    /// Return true if this Data object is a integer number.
    auto isInteger() const -> bool;

    /// Return true if this Data object is a boolean value.
    auto isBoolean() const -> bool;

    /// Return true if this Data object is a real object.
    auto isReal() const -> bool;

    /// Return true if this Data object is a Param object.
    auto isParam() const -> bool;

    /// Return true if this Data object is a dictionary object.
    auto isDict() const -> bool;

    /// Return true if this Data object is a list object.
    auto isList() const -> bool;

    /// Return true if this Data object is a null value.
    auto isNull() const -> bool;

    /// Return the child Data object with given key, presuming this Data object is a dictionary.
    auto operator[](Chars key) const -> Data const&;

    /// Return the child Data object with given key, presuming this Data object is a dictionary.
    auto operator[](String const& key) const -> Data const&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    auto operator[](int const& index) const -> Data const&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    auto operator[](Index const& index) const -> Data const&;

    /// Return the child Data object with given key, presuming this Data object is a dictionary.
    auto at(String const& key) const -> Data const&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    auto at(Index const& index) const -> Data const&;

    /// Return the child Data object with given key or create it if not existent, presuming this Data object is a dictionary.
    auto at(String const& key) -> Data&;

    /// Return the child Data object with given index or create it if not existent, presuming this Data object is a list.
    auto at(Index const& index) -> Data&;

    /// Return the child Data object whose `attribute` has a given `value`, presuming this Data object is a list.
    auto with(String const& attribute, String const& value) const -> Data const&;

    /// Add a Data object to this Data object, which becomes a list if not already.
    auto add(Data const& data) -> Data&;

    /// Add a Data object with given key to this Data object, which becomes a dictionary if not already.
    auto add(Chars key, Data const& data) -> Data&;

    /// Add a Data object with given key to this Data object, which becomes a dictionary if not already.
    auto add(String const& key, Data const& data) -> Data&;

    /// Return true if a child parameter exists with given key, presuming this Data object is a dictionary.
    auto exists(String const& key) const -> bool;

    /// Convert this Data object to a `bool` value; raises an error if not convertible to `bool`.
    operator bool() const;

    /// Convert this Data object to a String value; raises an error if not convertible to String.
    operator String() const;

    /// Convert this Data object to a Param value; raises an error if not convertible to Param.
    operator Param() const;

    /// Used to allow conversion of objects with custom types to Data objects.
    template<typename Type>
    struct Encode
    {
        /// Evaluate the conversion of an object with custom type to a Data object.
        static auto eval(Data& data, Type const& obj) -> void
        {
            errorif(true, "Cannot convert an object of type ", typeid(Type).name(), " to Data because Convert::encode was not defined for it.");
        }
    };

    /// Used to allow conversion of Data objects to objects with custom types.
    template<typename Type>
    struct Decode
    {
        /// Evaluate the conversion of a Data object to an object with custom type.
        static auto decode(Data const& data, Type& obj) -> void
        {
            errorif(true, "Cannot convert an object a Data object to an object of type ", typeid(Type).name(), " because Convert::decode was not defined for it.");
        }
    };

    template<typename T>
    auto operator=(T const& obj) -> Data&
    {
        Encode<T>::eval(*this, obj);
        return *this;
    }

    template<typename T>
    Data(T const& obj)
    {
        assign(obj);
    }

    template<typename T>
    auto assign(T const& obj) -> void
    {
        Encode<T>::eval(*this, obj);
    }

    auto assign(Chars obj) -> void { tree = String(obj); }

private:
    Any tree;
};

template<> inline auto Data::assign(bool const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(String const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(int const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(double const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(real const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(Param const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(Vec<Data> const& obj) -> void { tree = obj; }
template<> inline auto Data::assign(Map<String, Data> const& obj) -> void { tree = obj; }

#define REAKTORO_DATA_ENCODE_DECLARE(Type, ...)                                      \
    /** Used to encode/serialize an instance of type `Type` into a Data object. */   \
    template<__VA_ARGS__> struct Data::Encode<Type> { static auto eval(Data& data, const Type& obj) -> void; };

#define REAKTORO_DATA_ENCODE_DEFINE(Type, ...)                                       \
    /** Encode/serialize an instance of type `Type` into a Data object. */           \
    auto Data::Encode<Type>::eval(Data& data, const Type& obj) -> void

#define REAKTORO_DATA_DECODE_DECLARE(Type, ...)                                      \
    /** Used to decode/deserialize a Data object into an instance of type `Type`. */ \
    template<__VA_ARGS__> struct Data::Decode<Type> { static auto eval(const Data& data, Type& obj) -> void; };

#define REAKTORO_DATA_DECODE_DEFINE(Type, ...)                                       \
    /** Decode/deserialize a Data object into an instance of type `Type`. */         \
    auto Data::Decode<Type>::eval(const Data& data, Type& obj) -> void

#define REAKTORO_COMMA ,

REAKTORO_DATA_ENCODE_DECLARE(int);
REAKTORO_DATA_DECODE_DECLARE(int);

REAKTORO_DATA_ENCODE_DECLARE(Index);
REAKTORO_DATA_DECODE_DECLARE(Index);

REAKTORO_DATA_ENCODE_DECLARE(bool);
REAKTORO_DATA_DECODE_DECLARE(bool);

REAKTORO_DATA_ENCODE_DECLARE(Chars);
REAKTORO_DATA_DECODE_DECLARE(Chars);

REAKTORO_DATA_ENCODE_DECLARE(String);
REAKTORO_DATA_DECODE_DECLARE(String);

REAKTORO_DATA_ENCODE_DECLARE(Vec<T>, typename T);
REAKTORO_DATA_DECODE_DECLARE(Vec<T>, typename T);

REAKTORO_DATA_ENCODE_DECLARE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N);
REAKTORO_DATA_DECODE_DECLARE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N);

REAKTORO_DATA_ENCODE_DECLARE(Pair<A REAKTORO_COMMA B>, typename A, typename B);
REAKTORO_DATA_DECODE_DECLARE(Pair<A REAKTORO_COMMA B>, typename A, typename B);

REAKTORO_DATA_ENCODE_DECLARE(Map<K REAKTORO_COMMA T>, typename K, typename T);
REAKTORO_DATA_DECODE_DECLARE(Map<K REAKTORO_COMMA T>, typename K, typename T);

REAKTORO_DATA_ENCODE_DECLARE(Param);
REAKTORO_DATA_DECODE_DECLARE(Param);

REAKTORO_DATA_ENCODE_DECLARE(real);
REAKTORO_DATA_DECODE_DECLARE(real);

} // namespace Reaktoro
