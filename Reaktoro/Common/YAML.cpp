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

#include "YAML.hpp"

// C++ includes
#include <sstream>

namespace Reaktoro {

auto yaml::parse(const char* input) -> yaml
{
    return yaml(YAML::Load(input));
}

auto yaml::parse(const std::string& input) -> yaml
{
    return yaml(YAML::Load(input));
}

auto yaml::parse(std::istream& input) -> yaml
{
    return yaml(YAML::Load(input));
}

yaml::yaml()
: YAML::Node()
{}

yaml::yaml(const Node& node)
: YAML::Node(node)
{}

auto yaml::at(const std::string& key) const -> yaml
{
    auto child = (*this)[key];
    errorif(!child, "Could not get YAML child node with key `", key, "` in:\n", repr());
    return child;
}

auto yaml::repr() const -> std::string
{
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

} // namespace Reaktoro
