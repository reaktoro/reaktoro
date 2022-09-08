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

#include "Data.hpp"

// Third-party includes
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
using yaml = YAML::Node;
using json = nlohmann::json;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace {

/// Check if string `str` is a number.
/// @param str The string being checked
/// @param[out] result The number in `str` as a double value if it is indeed a number.
bool isNumber(String const& str, double& result)
{
    char* end;
    result = std::strtod(str.c_str(), &end);
    if(end == str.c_str() || *end != '\0')
        return false;
    return true;
}

auto convertYaml(yaml const& obj) -> Data;

auto convertYamlScalar(yaml const& obj) -> Data
{
    auto const& word = obj.Scalar();
    auto number = 0.0;
    if(isNumber(word, number))
        return Data(number);
    else
    {
        if(word == "true" || word == "True")
            return Data(true);
        if(word == "false" || word == "False")
            return Data(false);
        return Data(word);
    }
}

auto convertYamlSequence(yaml const& obj) -> Data
{
    Data result;
    for(auto i = 0; i < obj.size(); ++i)
        result.add(convertYaml(obj[i]));
    return result;
}

auto convertYamlMap(yaml const& obj) -> Data
{
    Data result;
    for(auto const& child : obj)
        result.add(child.first.as<String>(), convertYaml(child.second));
    return result;
}

auto convertYaml(yaml const& obj) -> Data
{
    switch(obj.Type())
    {
        case YAML::NodeType::Scalar: return convertYamlScalar(obj);
        case YAML::NodeType::Sequence: return convertYamlSequence(obj);
        case YAML::NodeType::Map: return convertYamlMap(obj);
        case YAML::NodeType::Null: return {};
    }

    errorif(true, "Could not convert YAML node to Data object: ", YAML::Dump(obj));

    return {};
}

auto convertJson(json const& obj) -> Data;

auto convertJsonObject(json const& obj) -> Data
{
    Data result;
    for(auto const& [key, value] : obj.items())
        result.add(key, convertJson(value));
    return result;
}

auto convertJsonArray(json const& obj) -> Data
{
    Data result;
    for(auto const& value : obj)
        result.add(convertJson(value));
    return result;
}

auto convertJson(json const& obj) -> Data
{
    switch(obj.type())
    {
        case json::value_t::null: return {};
        case json::value_t::object: return convertJsonObject(obj);
        case json::value_t::array: return convertJsonArray(obj);
        case json::value_t::string: return obj.get<String>();
        case json::value_t::boolean: return obj.get<bool>();
        case json::value_t::number_integer: return obj.get<int>();
        case json::value_t::number_unsigned: return obj.get<int>();
        case json::value_t::number_float: return obj.get<double>();
    }

    errorif(true, "Could not convert JSON node to Data object: ", obj.dump());

    return {};
}

/// An auxiliary type to change locale and ensure its return to original.
/// This is needed to avoid certain issues with pugixml related to how decimal numbers are represented in different languages.
struct ChangeLocale
{
    const String old_locale;

    explicit ChangeLocale(const char* new_locale) : old_locale(std::setlocale(LC_NUMERIC, nullptr))
    {
        std::setlocale(LC_NUMERIC, new_locale);
    }

    ~ChangeLocale()
    {
        std::setlocale(LC_NUMERIC, old_locale.c_str());
    }
};

} // namespace

Data::Data()
: tree(nullptr)
{
}

Data::Data(bool const& value)
: tree(value)
{
}

Data::Data(Chars const& value)
: tree(String(value))
{
}

Data::Data(String const& value)
: tree(value)
{
}

Data::Data(Param const& value)
: tree(value)
{
}

Data::Data(Map<String, Data> const& value)
: tree(value)
{
}

Data::Data(Vec<Data> const& value)
: tree(value)
{
}

Data::Data(yaml const& obj)
: Data(convertYaml(obj))
{
}

Data::Data(json const& obj)
: Data(convertJson(obj))
{
}

auto Data::fromYaml(Chars const& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return Data(YAML::Load(input));
}

auto Data::fromYaml(String const& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return Data(YAML::Load(input));
}

auto Data::fromYaml(std::istream& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return Data(YAML::Load(input));
}

auto Data::fromJson(Chars const& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return Data(nlohmann::json::parse(input));

}

auto Data::fromJson(String const& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return Data(nlohmann::json::parse(input));

}

auto Data::fromJson(std::istream& input) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    nlohmann::json obj;
    input >> obj;
    return Data(obj);
}

auto Data::string() const -> String const&
{
    return std::any_cast<String const&>(tree);
}

auto Data::number() const -> double
{
    return std::any_cast<Param const&>(tree).value();
}

auto Data::integer() const -> int
{
    return std::any_cast<Param const&>(tree).value();
}

auto Data::boolean() const -> bool
{
    return std::any_cast<bool const&>(tree);
}

auto Data::param() const -> Param const&
{
    return std::any_cast<Param const&>(tree);
}

auto Data::dict() const -> Map<String, Data> const&
{
    return std::any_cast<Map<String, Data> const&>(tree);
}

auto Data::list() const -> Vec<Data> const&
{
    return std::any_cast<Vec<Data> const&>(tree);
}

auto Data::null() const -> std::nullptr_t
{
    return std::any_cast<std::nullptr_t>(tree);
}

auto Data::isString() const -> bool
{
    return std::any_cast<String>(&tree);
}

auto Data::isNumber() const -> bool
{
    return std::any_cast<double>(&tree);
}

auto Data::isInteger() const -> bool
{
    return std::any_cast<int>(&tree);
}

auto Data::isBoolean() const -> bool
{
    return std::any_cast<bool>(&tree);
}

auto Data::isParam() const -> bool
{
    return std::any_cast<Param>(&tree);
}

auto Data::isDict() const -> bool
{
    return std::any_cast<Map<String, Data>>(&tree);
}

auto Data::isList() const -> bool
{
    return std::any_cast<Vec<Data>>(&tree);
}

auto Data::isNull() const -> bool
{
    return std::any_cast<std::nullptr_t>(&tree);
}

auto Data::operator[](String const& key) const -> Data const&
{
    auto const& obj = dict();
    auto const it = obj.find(key);
    errorif(it == obj.end(), "Could not find child data block with given key `,", key, "`.");
    return it->second;
}

auto Data::operator[](Index const& index) const -> Data const&
{
    auto const& vec = list();
    errorif(index >= vec.size(), "Could not retrieve child data block with index `,", index, "` because the list has size ", vec.size(), ".");
    return vec[index];
}

auto Data::with(String const& attribute, String const& value) const -> Data const&
{
    errorif(!isList(), "Expecting Data object to be a list when using Data::with method.");
    for(auto const& entry : list())
        if(entry[attribute].string() == value)
            return entry;
    errorif(true, "Could not find any data block whose attribute `", attribute, "` has value `", value, "`.");
    return {};
}

auto Data::add(Data const& data) -> Data&
{
    if(!isList())
        tree = Vec<Data>();
    auto& vec = std::any_cast<Vec<Data>&>(tree);
    vec.push_back(data);
    return vec.back();
}

auto Data::add(Chars const& key, Data const& data) -> Data&
{
    if(!isDict())
        tree = Map<String, Data>();
    auto& dict = std::any_cast<Map<String, Data>&>(tree);
    dict[key] = data;
    return dict[key];
}

auto Data::add(String const& key, Data const& data) -> Data&
{
    return add(key.c_str(), data);
}

auto Data::exists(String const& key) const -> bool
{
    if(!isDict())
        return false;
    auto obj = dict();
    return obj.find(key) != obj.end();
}

} // namespace Reaktoro

