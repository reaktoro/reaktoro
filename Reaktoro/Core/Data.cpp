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

// C++ includes
#include <fstream>

// Third-party includes
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
using yaml = YAML::Node;
using json = nlohmann::json;

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace {

// ==========================================================================================
// AUXILIARY METHODS
// ==========================================================================================

/// Check if string `str` is a number.
/// @param str The string being checked
/// @param[out] result The number in `str` as a double value if it is indeed a number.
bool isFloat(String const& str, double& result)
{
    auto const& newstr = lowercase(str);
    if(oneof(newstr, ".nan"))
    {
        result = std::numeric_limits<double>::quiet_NaN();
        return true;
    }
    if(oneof(newstr, "+inf", "+.inf", ".inf"))
    {
        result = std::numeric_limits<double>::infinity();
        return true;
    }
    if(oneof(newstr, "-inf", "-.inf"))
    {
        result = -std::numeric_limits<double>::infinity();
        return true;
    }

    if(oneof(newstr, "nan", "inf"))
        return false; // nan, inf, NaN, InF, etc. are considered words, not numbers! In fact, InF and NaF are valid substance names present in the NASA-CEA database!

    char* end;
    result = std::strtod(str.c_str(), &end);
    if(end == str.c_str() || *end != '\0')
        return false;
    return true;
}

// ==========================================================================================
// METHODS TO CONVERT YAML TO DATA
// ==========================================================================================

auto convertYamlToData(yaml const& obj) -> Data;

auto convertYamlScalarToData(yaml const& obj) -> Data
{
    auto const& word = obj.Scalar();
    auto number = 0.0;
    if(isFloat(word, number))
        return Param(number);
    else
    {
        if(word == "true" || word == "True")
            return true;
        if(word == "false" || word == "False")
            return false;
        return word;
    }
}

auto convertYamlSequenceToData(yaml const& obj) -> Data
{
    Data result;
    for(auto i = 0; i < obj.size(); ++i)
        result.add(convertYamlToData(obj[i]));
    return result;
}

auto convertYamlMapToData(yaml const& obj) -> Data
{
    Data result;
    for(auto const& child : obj)
        result.add(child.first.as<String>(), convertYamlToData(child.second));
    return result;
}

auto convertYamlToData(yaml const& obj) -> Data
{
    switch(obj.Type())
    {
        case YAML::NodeType::Scalar: return convertYamlScalarToData(obj);
        case YAML::NodeType::Sequence: return convertYamlSequenceToData(obj);
        case YAML::NodeType::Map: return convertYamlMapToData(obj);
        case YAML::NodeType::Null: return {};
    }

    errorif(true, "Could not convert YAML node to Data object: ", YAML::Dump(obj));

    return {};
}

// ==========================================================================================
// METHODS TO CONVERT JSON TO DATA
// ==========================================================================================

auto convertJsonToData(json const& obj) -> Data;

auto convertJsonObjectToData(json const& obj) -> Data
{
    Data result;
    for(auto const& [key, value] : obj.items())
        result.add(key, convertJsonToData(value));
    return result;
}

auto convertJsonArrayToData(json const& obj) -> Data
{
    Data result;
    for(auto const& value : obj)
        result.add(convertJsonToData(value));
    return result;
}

auto convertJsonToData(json const& obj) -> Data
{
    switch(obj.type())
    {
        case json::value_t::null: return {};
        case json::value_t::object: return convertJsonObjectToData(obj);
        case json::value_t::array: return convertJsonArrayToData(obj);
        case json::value_t::string: return obj.get<String>();
        case json::value_t::boolean: return obj.get<bool>();
        case json::value_t::number_integer: return obj.get<int>();
        case json::value_t::number_unsigned: return obj.get<int>();
        case json::value_t::number_float: return obj.get<double>();
    }

    errorif(true, "Could not convert JSON node to Data object: ", obj.dump());

    return {};
}

// ==========================================================================================
// METHODS TO CONVERT DATA TO YAML AND JSON
// ==========================================================================================

template<typename Format>
auto convertDataTo(Data const& data) -> Format;

template<typename Format>
auto convertDataDictTo(Data const& data) -> Format
{
    assert(data.isDict());
    Format res;
    for(auto const& [key, value] : data.asDict())
        res[key] = convertDataTo<Format>(value);
    return res;
}

template<typename Format>
auto convertDataListTo(Data const& data) -> Format
{
    assert(data.isList());
    Format res;
    for(auto const& value : data.asList())
        res.push_back(convertDataTo<Format>(value));
    return res;
}

template<typename Format>
auto convertDataTo(Data const& data) -> Format
{
    if(data.isNull()) return Format();
    if(data.isBoolean()) return Format(data.asBoolean());
    if(data.isString()) return Format(data.asString());
    if(data.isParam()) return Format(data.asFloat());
    if(data.isDict()) return convertDataDictTo<Format>(data);
    if(data.isList()) return convertDataListTo<Format>(data);
    errorif(true, "Could not convert this Data object to an YAML or JSON as the Data object is not in a valid state.");
    return {};
}

auto convertDataToYaml(Data const& data) -> yaml
{
    return convertDataTo<yaml>(data);
}

auto convertDataToJson(Data const& data) -> json
{
    return convertDataTo<json>(data);
}

// ==========================================================================================
// CLASS TO ENSURE A COMMON LOCALE IS KEPT WHEN DEALING WITH YAML AND JSON
// ==========================================================================================

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

// ==========================================================================================
// IMPLEMENTATION OF CLASS DATA
// ==========================================================================================

Data::Data()
: tree(nullptr)
{
}

auto Data::parse(Chars text) -> Data
{
    return parseYaml(text);
}

auto Data::parse(String const& text) -> Data
{
    return parseYaml(text);
}

auto Data::parse(std::istream& text) -> Data
{
    return parseYaml(text);
}

auto Data::parseYaml(Chars text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return convertYamlToData(YAML::Load(text));
}

auto Data::parseYaml(String const& text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return convertYamlToData(YAML::Load(text));
}

auto Data::parseYaml(std::istream& text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return convertYamlToData(YAML::Load(text));
}

auto Data::parseJson(Chars text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return convertJsonToData(nlohmann::json::parse(text));
}

auto Data::parseJson(String const& text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return convertJsonToData(nlohmann::json::parse(text));
}

auto Data::parseJson(std::istream& text) -> Data
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    nlohmann::json obj;
    text >> obj;
    return convertJsonToData(obj);
}

auto Data::load(String const& path) -> Data
{
    const auto words = split(path, ".");
    if(words.back() == "yaml" || words.back() == "yml")
        return Data::loadYaml(path);
    if(words.back() == "json")
        return Data::loadJson(path);
    errorif(true, "The given path `", path, "` was expected to point to a file containing one of the following file extensions: yml, yaml, json. Please rename your file accordingly.");
    return {};
}

auto Data::loadYaml(String const& path) -> Data
{
    std::ifstream f(path);
    errorif(f.fail(), "There was an error finding your YAML file at `", path, "`. Ensure this file exists and prefer global path strings such as \"/home/mary/data.json\" in Linux and macOS or \"C:\\\\Users\\\\Mary\\\\data.json\" in Windows.");
    yaml doc;
    try { doc = YAML::Load(f); }
    catch(std::exception e)
        errorif(true, "There was an error parsing your YAML file at `", path, "`. Ensure this file is properly formatted (e.g., inconsistent indentation). Try using some online YAML validator to find the error. More details about the error below:\n\n", e.what());
    return convertYamlToData(doc);
}

auto Data::loadJson(String const& path) -> Data
{
    std::ifstream f(path);
    errorif(f.fail(), "There was an error finding your JSON file at `", path, "`. Ensure this file exists and prefer global path strings such as \"/home/mary/data.json\" in Linux and macOS or \"C:\\\\Users\\\\Mary\\\\data.json\" in Windows.");
    json doc;
    try { doc = json::parse(f); }
    catch(std::exception e)
        errorif(true, "There was an error parsing your JSON file at `", path, "`. Ensure this file is properly formatted (e.g., missing closing brackets). Try using some online JSON validator to find the error. More details about the error below:\n\n", e.what());
    return convertJsonToData(doc);
}


auto Data::asString() const -> String const&
{
    errorif(!isString(), "Cannot convert this Data object to a String.");
    return std::any_cast<String const&>(tree);
}

auto Data::asBoolean() const -> bool
{
    errorif(!isBoolean(), "Cannot convert this Data object to a boolean value.");
    return std::any_cast<bool const&>(tree);
}

auto Data::asInteger() const -> int
{
    if(isParam())
        return std::any_cast<Param const&>(tree).value().val();
    else errorif(true, "Cannot convert this Data object to an integer number. This Data object should be a Param object.");
}

auto Data::asFloat() const -> double
{
    if(isParam())
        return std::any_cast<Param const&>(tree).value().val();
    else errorif(true, "Cannot convert this Data object to a float number. This Data object should be a Param object.");
}

auto Data::asReal() const -> real const&
{
    if(isParam())
        return std::any_cast<Param const&>(tree).value();
    else errorif(true, "Cannot convert this Data object to a real number. This Data object should be a Param object.");
}

auto Data::asParam() const -> Param const&
{
    errorif(!isParam(), "Cannot convert this Data object to a Param object. This Data object should be a Param object.");
    return std::any_cast<Param const&>(tree);
}

auto Data::asDict() const -> Dict<String, Data> const&
{
    errorif(!isDict(), "Cannot convert this Data object to a dictionary object.");
    return std::any_cast<Dict<String, Data> const&>(tree);
}

auto Data::asList() const -> Vec<Data> const&
{
    errorif(!isList(), "Cannot convert this Data object to a list object.");
    return std::any_cast<Vec<Data> const&>(tree);
}

auto Data::isBoolean() const -> bool
{
    return std::any_cast<bool>(&tree);
}

auto Data::isString() const -> bool
{
    return std::any_cast<String>(&tree);
}

auto Data::isParam() const -> bool
{
    return std::any_cast<Param>(&tree);
}

auto Data::isDict() const -> bool
{
    return std::any_cast<Dict<String, Data>>(&tree);
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
    return at(key);
}

auto Data::operator[](Index const& index) const -> Data const&
{
    return at(index);
}

auto Data::operator[](String const& key) -> Data&
{
    if(isNull())
        tree = Dict<String, Data>();
    errorif(!isDict(), "Methods Data::at(key) and Data::operator[key], with key `", key, "` can only be used when the Data object is a dictionary.");
    auto& obj = std::any_cast<Dict<String, Data>&>(tree);
    return obj[key];
}

auto Data::operator[](Index const& index) -> Data&
{
    if(isNull() && index == 0)
        tree = Vec<Data>();
    errorif(!isList(), "Methods Data::at(index) and Data::operator[index] can only be used when the Data object is a list.");
    auto& list = std::any_cast<Vec<Data>&>(tree);
    errorif(index >= list.size(), "Could not retrieve data block with index ", index, " because the list has size ", list.size(), ".");
    return list[index];
}

auto Data::at(String const& key) const -> Data const&
{
    errorif(!isDict(), "Methods Data::at(key) and Data::operator[key], with key `", key, "` can only be used when the Data object (const in this context) is a dictionary.");
    auto const& dict = std::any_cast<Dict<String, Data> const&>(tree);
    auto const it = dict.find(key);
    errorif(it == dict.end(), "Could not find data block with given key `", key, "`.");
    return it->second;
}

auto Data::at(Index const& index) const -> Data const&
{
    errorif(!isList(), "Methods Data::at(index) and Data::operator[index] can only be used when the Data object (const in this context) is a list.");
    auto const& list = std::any_cast<Vec<Data> const&>(tree);
    errorif(index >= list.size(), "Could not retrieve data block with index ", index, " because the list has size ", list.size(), ".");
    return list[index];
}

auto Data::optional(String const& key) const -> Opt
{
    errorif(!isDict(), "Method Data::optional(key), with key `", key, "` can only be used when the Data object is a dictionary.");
    auto const& dict = std::any_cast<Dict<String, Data> const&>(tree);
    auto const it = dict.find(key);
    return it != dict.end() ? Opt{&it->second} : Opt{};
}

auto Data::required(String const& key) const -> Data const&
{
    errorif(!isDict(), "Method Data::required(key), with key `", key, "` can only be used when the Data object is a dictionary.");
    auto const& dict = std::any_cast<Dict<String, Data> const&>(tree);
    auto const it = dict.find(key);
    errorif(it == dict.end(), "Could not find required data block with key `", key, "` in the Data object.");
    return it->second;
}

auto Data::with(String const& attribute, String const& value) const -> Data const&
{
    errorif(!isList(), "Expecting Data object to be a list when using Data::with method.");
    for(auto const& entry : asList())
        if(entry[attribute].asString() == value)
            return entry;
    errorif(true, "Could not find any data block whose attribute `", attribute, "` has value `", value, "`.");
    return *this;
}

auto Data::add(Data const& value) -> Data&
{
    if(isNull())
        tree = Vec<Data>();
    errorif(!isList(), "Method Data::add(value) can only be used when the Data object is a list or null.");
    auto& list = std::any_cast<Vec<Data>&>(tree);
    list.push_back(value);
    return list.back();
}

auto Data::add(Chars key, Data const& data) -> Data&
{
    if(isNull())
        tree = Dict<String, Data>();
    errorif(!isDict(), "Method Data::add(key, value) can only be used when the Data object is a dictionary or null.");
    auto& dict = std::any_cast<Dict<String, Data>&>(tree);
    dict[key] = data;
    return dict[key];
}

auto Data::add(String const& key, Data const& data) -> Data&
{
    return add(key.c_str(), data);
}

auto Data::reset() -> void
{
    tree = nullptr;
}

auto Data::exists(String const& key) const -> bool
{
    if(!isDict())
        return false;
    auto obj = asDict();
    return obj.find(key) != obj.end();
}

auto Data::dumpYaml() const -> String
{
    yaml doc = convertDataToYaml(*this);
    return YAML::Dump(doc);
}

auto Data::dumpJson() const -> String
{
    json doc = convertDataToJson(*this);
    return doc.dump();
}

auto Data::repr() const -> String
{
    return dumpYaml();
}

} // namespace Reaktoro

