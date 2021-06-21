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
namespace {

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

} // namespace (anonymous)

auto yaml::parse(const char* input) -> yaml
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return yaml(YAML::Load(input));
}

auto yaml::parse(const std::string& input) -> yaml
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
    return yaml(YAML::Load(input));
}

auto yaml::parse(std::istream& input) -> yaml
{
    const auto guard = ChangeLocale("C"); // Change locale to C before parsing (this is reset at destruction of `guard`).
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
    if(IsNull()) ss << "null";
    else ss << *this;
    return ss.str();
}

} // namespace Reaktoro
