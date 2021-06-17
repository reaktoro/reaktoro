// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

// C++ includes
#include <istream>
#include <optional>
#include <string>

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

/// A class used to serialize/deserialize other types based on yaml format.
class yaml : public YAML::Node
{
public:
    /// Parse a given input string in YAML format.
    static auto parse(const char* input) -> yaml;

    /// Parse a given input string in YAML format.
    static auto parse(const std::string& input) -> yaml;

    /// Parse a given input stream in YAML format.
    static auto parse(std::istream& input) -> yaml;

    /// Construct an object of type yaml.
    yaml();

    /// Construct an object of type yaml from a given Node object.
    yaml(const YAML::Node& node);

    /// Construct an object of type yaml from a given object of type @p T.
    template<typename T>
    yaml(const T& obj) : YAML::Node(obj) {}

    /// Append a child node with a given value only if value is not default value.
    template<typename T, typename U = T>
    auto appendIfNotDefault(const std::string& key, const T& value, const U& defaultval = U{}) -> void {
        if(value != defaultval) (*this)[key] = value; }

    /// Return a child node with given key if found, otherwise raise an error.
    auto at(const std::string& key) const -> yaml;

    /// Return the value of the node if not empty, otherwise a given fallback value.
    template<typename T>
    auto value(const T& fallback = T{}) const -> T { return IsDefined() ? as<T>() : fallback; }

    /// Return a string representation of this yaml object.
    auto repr() const -> std::string;

    /// Return child with given key.
    template<typename Key>
    auto operator[](const Key& key) const -> const yaml { return YAML::Node::operator[](key); }

    /// Return child with given key.
    template<typename Key>
    auto operator[](const Key& key) -> yaml { return YAML::Node::operator[](key); }

    /// Assign this yaml node with given value.
    template<typename T>
    auto operator=(const T& value) { YAML::Node::operator=(value); return *this; }

    /// Implicitly convert this yaml object into another type.
    template<typename T>
    operator T() const { return as<T>(); }

    /// Transfer the value at this yaml node to argument @p value.
    template<typename T>
    auto to(T& value) const -> void
    {
        try {
            value = as<T>();
        }
        catch(...) {
            errorif(true, "Could not convert YAML node to requested value type. The node is:\n", repr());
        }
    }

    /// Transfer the value at this yaml node to argument @p value. Use fallback in case there of empty node.
    template<typename T, typename U>
    auto to(T& value, const U& fallback) const -> void
    {
        if(!IsDefined()) value = fallback;
        else to(value);
    }

    /// Transfer the value at optional child node with @p key to optional object @p obj.
    template<typename T>
    auto copyOptionalChildValueTo(const std::string& key, T& obj, const T& defaultval) const
    {
        auto child = (*this)[key];
        if(child.IsDefined())
            child.to<T>(obj);
        else obj = defaultval;
    }

    /// Transfer the value at required child node with @p key to object @p obj.
    template<typename T>
    auto copyRequiredChildValueTo(const std::string& key, T& obj) const
    {
        auto child = (*this)[key];
        errorif(!child.IsDefined(), "Could not find a required YAML node with key `", key, "` in parent node:\n", repr());
        child.to(obj);
    }

    template<typename T>
    struct encode;

    template<typename T>
    struct decode;
};

#define REAKTORO_YAML_ENCODE_DECLARE(Type) \
    template<> struct yaml::encode<Type> { static auto eval(yaml& node, const Type& obj) -> void; };

#define REAKTORO_YAML_ENCODE_DEFINE(Type) \
    auto yaml::encode<Type>::eval(yaml& node, const Type& obj) -> void

#define REAKTORO_YAML_DECODE_DECLARE(Type) \
    template<> struct yaml::decode<Type> { static auto eval(const yaml& node, Type& obj) -> void; };

#define REAKTORO_YAML_DECODE_DEFINE(Type) \
    auto yaml::decode<Type>::eval(const yaml& node, Type& obj) -> void

} // namespace Reaktoro

namespace YAML {

template<typename Type>
struct convert
{
    static auto encode(const Type& obj)
    {
        Reaktoro::yaml node;
        Reaktoro::yaml::encode<Type>::eval(node, obj);
        return node;
    }

    static auto decode(const Node& node, Type& obj)
    {
        Reaktoro::yaml::decode<Type>::eval(node, obj);
        return true;
    }
};

} // namespace YAML
