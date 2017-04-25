// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <unsupported/cpp-interpreter/Yaml.hpp>

namespace Reaktoro {

auto str(const Node& node) -> std::string
{
    std::stringstream ss; ss << node;
    return ss.str();
}

auto keynode(const Node& node) -> Node
{
    Assert(node.IsMap(), "Could not get the key of node `" + str(node) + "`.",
        "This is only allowed for nodes that are maps.");
    return node.begin()->first;
}

auto valnode(const Node& node) -> Node
{
    Assert(node.IsMap(), "Could not get the value of node `" + str(node) + "`.",
        "This is only allowed for nodes that are maps.");
    return node.begin()->second;
}

auto keyword(const Node& node) -> std::string
{
    std::string key = str(keynode(node));
    auto words = split(key);
    return words[0];
}

auto identifier(const Node& node) -> std::string
{
    std::string key = str(keynode(node));
    auto words = split(key);
    return words.size() > 1 ? words[1] : "";
}

auto operator+(std::string str, const Node& node) -> std::string
{
    std::stringstream res;
    res << str << node;
    return res.str();
}

auto operator+(const Node& node, std::string str) -> std::string
{
    std::stringstream res;
    res << node << str;
    return res.str();
}

} // namespace Reaktoro
