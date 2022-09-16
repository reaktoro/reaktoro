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
#include <tsl/ordered_map.h>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Core/Param.hpp>

namespace Reaktoro {

template<class Key, class T>
using Dict = tsl::ordered_map<Key, T>;

} // namespace Reaktoro

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

    /// Return a Data object by parsing an YAML formatted string.
    static auto parse(Chars text) -> Data;

    /// Return a Data object by parsing an YAML formatted string.
    static auto parse(String const& text) -> Data;

    /// Return a Data object by parsing an YAML formatted string.
    static auto parse(std::istream& text) -> Data;

    /// Return a Data object by parsing an YAML formatted string.
    static auto parseYaml(Chars text) -> Data;

    /// Return a Data object by parsing an YAML formatted string.
    static auto parseYaml(String const& text) -> Data;

    /// Return a Data object by parsing an YAML formatted string.
    static auto parseYaml(std::istream& text) -> Data;

    /// Return a Data object by parsing a JSON formatted string.
    static auto parseJson(Chars text) -> Data;

    /// Return a Data object by parsing a JSON formatted string.
    static auto parseJson(String const& text) -> Data;

    /// Return a Data object by parsing a JSON formatted string.
    static auto parseJson(std::istream& text) -> Data;

    /// Return a Data object by parsing either an YAML or JSON formatted file at a given path.
    /// Ensure `path` terminates with extensions `.yml`, `.yaml`, or `.json`, which are used to identify the file format.
    static auto load(String const& path) -> Data;

    /// Return a Data object by parsing an YAML formatted file at a given path.
    static auto loadYaml(String const& path) -> Data;

    /// Return a Data object by parsing a JSON formatted file at a given path.
    static auto loadJson(String const& path) -> Data;

    /// Return this Data object as a boolean value.
    auto asBoolean() const -> bool;

    /// Return this Data object as a string.
    auto asString() const -> String const&;

    /// Return this Data object as an integer number.
    auto asInteger() const -> int;

    /// Return this Data object as a float number.
    auto asFloat() const -> double;

    /// Return this Data object as a Param object.
    auto asParam() const -> Param const&;

    /// Return this Data object as a dictionary object.
    auto asDict() const -> Dict<String, Data> const&;

    /// Return this Data object as a list object.
    auto asList() const -> Vec<Data> const&;

    /// Return this Data object as a nullptr value.
    auto asNull() const -> Nullptr;

    /// Return true if this Data object is a boolean value.
    auto isBoolean() const -> bool;

    /// Return true if this Data object is a string.
    auto isString() const -> bool;

    /// Return true if this Data object is an integer number.
    auto isInteger() const -> bool;

    /// Return true if this Data object is a float number.
    auto isFloat() const -> bool;

    /// Return true if this Data object is a Param object.
    auto isParam() const -> bool;

    /// Return true if this Data object is a dictionary object.
    auto isDict() const -> bool;

    /// Return true if this Data object is a list object.
    auto isList() const -> bool;

    /// Return true if this Data object is a null value.
    auto isNull() const -> bool;

    /// Return the child Data object with given key, presuming this Data object is a dictionary.
    /// This method throws an error if this Data object is not a dictionary or if the given key does not exist.
    auto operator[](String const& key) const -> Data const&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    /// This method throws an error if this Data object is not a list or if the given index is out of bounds.
    auto operator[](Index const& index) const -> Data const&;

    /// Return the child Data object with given key if, presuming this Data object is a dictionary.
    /// If this Data object is null, this method converts it into a dictionary with one entry with given key whose associated value is null.
    /// If this Data object is a dictionary without an entry with given key, a new entry with given key is created whose associated value is null.
    /// Otherwise, a runtime error is thrown.
    auto operator[](String const& key) -> Data&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    /// If this Data object is null and `index` is 0, this method converts it into a list with one entry whose value is null.
    /// Otherwise, this method throws an error if this Data object is not a list or it is a list and given index is out of bounds.
    auto operator[](Index const& index) -> Data&;

    /// Return the child Data object with given key, presuming this Data object is a dictionary.
    /// This method throws an error if this Data object is not a dictionary or if the given key does not exist.
    auto at(String const& key) const -> Data const&;

    /// Return the child Data object with given index, presuming this Data object is a list.
    /// This method throws an error if this Data object is not a list or if the given index is out of bounds.
    auto at(Index const& index) const -> Data const&;

    /// Used as the return type of method Data::optional.
    struct Opt
    {
        Data const*const ptrdata = nullptr;

        /// Decode this Opt object into an object of type `T`, or leave it unchanged if no data is available.
        template<typename T>
        auto to(T& obj) const -> void
        {
            if(ptrdata)
                obj = ptrdata->as<T>();
        }
    };

    /// Return an optional child Data object with given key, presuming this Data object is a dictionary.
    /// This method throws an error if this Data object is not a dictionary.
    auto optional(String const& key) const -> Opt;

    /// Return a required to exist child Data object with given key, presuming this Data object is a dictionary.
    /// This method throws an error if this Data object is not a dictionary or if the given key does not exist.
    auto required(String const& key) const -> Data const&;

    /// Return the child Data object whose `attribute` has a given `value`, presuming this Data object is a list.
    auto with(String const& attribute, String const& value) const -> Data const&;

    /// Add a Data object to this Data object, which becomes a list if not already.
    auto add(Data data) -> void;

    /// Add a Data object with given key to this Data object, which becomes a dictionary if not already.
    auto add(String const& key, Data data) -> void;

    /// Update this Data object with data given in another, with key-value pairs being either added or overwritten.
    auto update(Data const& data) -> void;

    /// Reset this Data object to a null state, deleting its current stored data.
    auto reset() -> void;

    /// Return true if a child Data object exists with given key, presuming this Data object is a dictionary.
    auto exists(String const& key) const -> bool;

    /// Return a YAML formatted string representing the state of this Data object.
    auto dump() const -> String;

    /// Return a YAML formatted string representing the state of this Data object.
    auto dumpYaml() const -> String;

    /// Return a JSON formatted string representing the state of this Data object.
    auto dumpJson() const -> String;

    /// Save the state of this Data object into a YAML formatted file.
    auto save(String const& filepath) const -> void;

    /// Save the state of this Data object into a YAML formatted file.
    auto saveYaml(String const& filepath) const -> void;

    /// Save the state of this Data object into a JSON formatted file.
    auto saveJson(String const& filepath) const -> void;

    /// Return a YAML formatted string representing the state of this Data object.
    auto repr() const -> String;

    /// Used to allow conversion of objects with custom types to Data objects.
    template<typename T>
    struct Encode
    {
        /// Evaluate the conversion of an object with custom type to a Data object.
        static auto eval(Data& data, T const& obj) -> void
        {
            errorif(true, "Cannot convert an object of type ", typeid(T).name(), " to Data because Encode::eval was not defined for it.");
        }
    };

    /// Used to allow conversion of Data objects to objects with custom types.
    template<typename T>
    struct Decode
    {
        /// Evaluate the conversion of a Data object to an object with custom type.
        static auto eval(Data const& data, T& obj) -> void
        {
            errorif(true, "Cannot convert an object a Data object to an object of type ", typeid(T).name(), " because Decode::eval was not defined for it.");
        }
    };

    /// Assign an object of type `T` to this Data object.
    template<typename T>
    auto operator=(T const& obj) -> Data&
    {
        assign(obj);
        return *this;
    }

    /// Convert this Data object to one of type `T`.
    template<typename T>
    explicit operator T() const
    {
        return as<T>();
    }

    /// Construct a Data object from one of type `T`.
    template<typename T>
    Data(T const& obj)
    {
        assign(obj);
    }

    /// Assign an object of type `T` to this Data object.
    template<typename T>
    auto assign(T const& obj) -> void
    {
        if constexpr(isOneOf<T, bool, int, double, String, Param, Vec<Data>, Dict<String, Data>, Nullptr>)
            tree = obj;
        else if constexpr(Reaktoro::isInteger<T>)
            tree = static_cast<int>(obj);
        else if constexpr(isFloatingPoint<T> || isSame<T, real>)
            tree = static_cast<double>(obj);
        else {
            reset();
            Encode<T>::eval(*this, obj);
        }
    }

    /// Assign a char to this Data object.
    auto assign(char obj) -> void
    {
        tree = String(1, obj);
    }

    /// Assign a raw string to this Data object.
    auto assign(Chars obj) -> void
    {
        tree = String(obj);
    }

    /// Convert this Data object into an object of type `T`.
    template<typename T>
    auto as() const -> T
    {
        if constexpr(isOneOf<T, bool, String, Param, Vec<Data>, Dict<String, Data>>) {
            const bool convertable = std::any_cast<T>(&tree);
            errorif(!convertable, "Could not convert from Data object to an object of type ", typeid(T).name(), " because this is not the type of the data stored nor it is convertible to that type.");
            return std::any_cast<T const&>(tree);
        }
        if constexpr(Reaktoro::isInteger<T>)
            return asInteger();
        if constexpr(isFloatingPoint<T> || isSame<T, real>)
            return asFloat();
        else {
            T obj;
            Decode<T>::eval(*this, obj);
            return obj;
        }
    }

    /// Decode this Data object into an object of type `T`.
    template<typename T>
    auto to(T& obj) const -> void
    {
        obj = as<T>();
    }

private:
    Any tree;
};

#define REAKTORO_DATA_ENCODE_DECLARE(T, ...)                                      \
    /** Used to encode/serialize an instance of type `T` into a Data object. */   \
    template<__VA_ARGS__> struct Data::Encode<T> { static auto eval(Data& data, const T& obj) -> void; };

#define REAKTORO_DATA_ENCODE_DEFINE(T, ...)                                       \
    /** Encode/serialize an instance of type `T` into a Data object. */           \
    auto Data::Encode<T>::eval(Data& data, const T& obj) -> void

#define REAKTORO_DATA_DECODE_DECLARE(T, ...)                                      \
    /** Used to decode/deserialize a Data object into an instance of type `T`. */ \
    template<__VA_ARGS__> struct Data::Decode<T> { static auto eval(const Data& data, T& obj) -> void; };

#define REAKTORO_DATA_DECODE_DEFINE(T, ...)                                       \
    /** Decode/deserialize a Data object into an instance of type `T`. */         \
    auto Data::Decode<T>::eval(const Data& data, T& obj) -> void

#define REAKTORO_COMMA ,

REAKTORO_DATA_ENCODE_DECLARE(Vec<T>, typename T);
REAKTORO_DATA_DECODE_DECLARE(Vec<T>, typename T);

REAKTORO_DATA_ENCODE_DECLARE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N);
REAKTORO_DATA_DECODE_DECLARE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N);

REAKTORO_DATA_ENCODE_DECLARE(Pair<A REAKTORO_COMMA B>, typename A, typename B);
REAKTORO_DATA_DECODE_DECLARE(Pair<A REAKTORO_COMMA B>, typename A, typename B);

REAKTORO_DATA_ENCODE_DECLARE(Map<K REAKTORO_COMMA T>, typename K, typename T);
REAKTORO_DATA_DECODE_DECLARE(Map<K REAKTORO_COMMA T>, typename K, typename T);

REAKTORO_DATA_ENCODE_DECLARE(Dict<K REAKTORO_COMMA T>, typename K, typename T);
REAKTORO_DATA_DECODE_DECLARE(Dict<K REAKTORO_COMMA T>, typename K, typename T);

template<typename T>
REAKTORO_DATA_ENCODE_DEFINE(Vec<T>, typename T)
{
    for(auto const& x : obj)
        data.add(x);
}

template<typename T>
REAKTORO_DATA_DECODE_DEFINE(Vec<T>, typename T)
{
    for(auto const& x : data.asList())
        obj.push_back(x.as<T>());
}

template<typename T, std::size_t N>
REAKTORO_DATA_ENCODE_DEFINE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N)
{
    for(auto const& x : obj)
        data.add(x);
}

template<typename T, std::size_t N>
REAKTORO_DATA_DECODE_DEFINE(Array<T REAKTORO_COMMA N>, typename T, std::size_t N)
{
    auto i = 0;
    for(auto const& x : data.asList())
        obj[i++] = x.as<T>();
}

template<typename A, typename B>
REAKTORO_DATA_ENCODE_DEFINE(Pair<A REAKTORO_COMMA B>, typename A, typename B)
{
    data.add(obj.first);
    data.add(obj.second);
}

template<typename A, typename B>
REAKTORO_DATA_DECODE_DEFINE(Pair<A REAKTORO_COMMA B>, typename A, typename B)
{
    auto const& l = data.asList();
    errorif(l.size() != 2, "Converting from Data to Pair requires the Data object to be a list with two entries.");
    obj.first = l[0].as<A>();
    obj.second = l[1].as<B>();
}

template<typename K, typename T>
REAKTORO_DATA_ENCODE_DEFINE(Map<K REAKTORO_COMMA T>, typename K, typename T)
{
    for(auto const& [k, v] : obj)
        data.add(k, v);
}

template<typename K, typename T>
REAKTORO_DATA_DECODE_DEFINE(Map<K REAKTORO_COMMA T>, typename K, typename T)
{
    for(auto const& [k, v] : data.asDict())
        obj[k] = v.template as<T>();
}

template<typename K, typename T>
REAKTORO_DATA_ENCODE_DEFINE(Dict<K REAKTORO_COMMA T>, typename K, typename T)
{
    for(auto const& [k, v] : obj)
        data.add(k, v);
}

template<typename K, typename T>
REAKTORO_DATA_DECODE_DEFINE(Dict<K REAKTORO_COMMA T>, typename K, typename T)
{
    for(auto const& [k, v] : data.asDict())
        obj[k] = v.template as<T>();
}

} // namespace Reaktoro
