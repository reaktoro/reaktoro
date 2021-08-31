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

#include "Params.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

Params::Params()
{}

auto Params::size() const -> Index
{
    auto count = 0;
    for(auto const& [key, val] : tree)
        if(val.type() == typeid(Params))
            count += std::any_cast<const Params&>(val).size();
        else count += 1;
    return count;
}

auto Params::at(const String& key) const -> const Params&
{
    const auto it = tree.find(key);
    error(it == tree.end(), "Could not find a node in the Params object with key `", key, "`.");
    return std::any_cast<const Params&>(it->second);
}

auto Params::get(const String& key) const -> const Param&
{
    error(tree.empty(), "Could not find a parameter in the empty Params object with key `", key, "`.");
    const auto it = tree.find(key);
    error(it == tree.end(), "Could not find a parameter in the Params object with key `", key, "`.");
    return std::any_cast<const Param&>(it->second);
}

auto Params::exists(const String& key) const -> bool
{
    const auto it = tree.find(key);
    return it != tree.end();
}

auto Params::set(const String& key, const Params& node) -> void
{
    tree[key] = node;
}

auto Params::set(const String& key, const Param& param) -> void
{
    tree[key] = param;
}

} // namespace Reaktoro

