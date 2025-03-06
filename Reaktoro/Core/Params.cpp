// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "Params.hpp"

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Embedded.hpp>

namespace Reaktoro {
namespace {

/// Create a Data object either from YAML or JSON text.
template<typename Text>
auto createDataFromYamlOrJson(String const& path, Text&& text) -> Data
{
    const auto words = split(path, ".");
    if(words.back() == "yaml" || words.back() == "yml")
        return Data::parseYaml(text);
    if(words.back() == "json")
        return Data::parseJson(text);
    errorif(true, "The given file path ", path, " does not indicate whether the file is YAML or JSON (expecting yaml, yml, or json as file extensions).");
    return {};
}

} // namespace

Params::Params()
{}

Params::Params(Data const& params)
{
    append(params);
}

auto Params::fromFile(String const& path) -> Params
{
    std::ifstream file(path);
    errorif(!file, "Could not load local parameter file with given path: ", path);
    return Params(createDataFromYamlOrJson(path, file));
}

auto Params::fromEmbeddedFile(String const& path) -> Params
{
    const String text = Embedded::get("params/" + path);
    return Params(createDataFromYamlOrJson(path, text));
}

auto Params::fromStringYAML(String const& text) -> Params
{
    auto data = Data::parseYaml(text);
    return Params(data);
}

auto Params::fromStringJSON(String const& text) -> Params
{
    auto data = Data::parseJson(text);
    return Params(data);
}

auto Params::embedded(String const& path) -> Params
{
    return fromEmbeddedFile(path);
}

auto Params::local(String const& path) -> Params
{
    return fromFile(path);
}

auto Params::data() const -> Data const&
{
    return m_data;
}

auto Params::append(Params const& other) -> Params&
{
    return append(other.data());
}

auto Params::append(Data const& other) -> Params&
{
    errorif(!other.isDict(), "Expecting Data object of dictionary type in Params object.");
    m_data.update(other);
    return *this;
}

auto Params::operator+=(Params const& other) -> Params&
{
    return append(other);
}

auto Params::operator[](String const& name) const -> Data const&
{
    return m_data[name];
}

} // namespace Reaktoro
