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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <any>
#include <cstddef>

namespace Reaktoro {

Data::Data()
: tree(nullptr)
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

Data::Data(double const& value)
: tree(value)
{
}

Data::Data(int const& value)
: tree(value)
{
}

Data::Data(bool const& value)
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

auto Data::string() const -> String const&
{
    return std::any_cast<String const&>(tree);
}

auto Data::number() const -> double const&
{
    return std::any_cast<double const&>(tree);
}

auto Data::integer() const -> int const&
{
    return std::any_cast<int const&>(tree);
}

auto Data::boolean() const -> bool const&
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

// Data::Data()
// {}

// auto Data::size() const -> Index
// {
//     auto count = 0;
//     for(auto const& [key, val] : tree)
//         if(val.type() == typeid(Data))
//             count += std::any_cast<const Data&>(val).size();
//         else count += 1;
//     return count;
// }

// auto Data::at(const String& key) const -> const Data&
// {
//     const auto it = tree.find(key);
//     error(it == tree.end(), "Could not find a node in the Data dict with key `", key, "`.");
//     return std::any_cast<const Data&>(it->second);
// }

// auto Data::get(const String& key) const -> const Param&
// {
//     error(tree.empty(), "Could not find a parameter in the empty Data dict with key `", key, "`.");
//     const auto it = tree.find(key);
//     error(it == tree.end(), "Could not find a parameter in the Data dict with key `", key, "`.");
//     return std::any_cast<const Param&>(it->second);
// }

// auto Data::exists(const String& key) const -> bool
// {
//     const auto it = tree.find(key);
//     return it != tree.end();
// }

// auto Data::set(const String& key, const Data& node) -> void
// {
//     tree[key] = node;
// }

// auto Data::set(const String& key, const Param& param) -> void
// {
//     tree[key] = param;
// }

} // namespace Reaktoro

