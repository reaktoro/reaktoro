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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Param.hpp>

namespace Reaktoro {

/// The class used to store and retrieve thermodynamic parameters.
/// @ingroup Core
class Params
{
public:
    /// Construct a default Params instance.
    Params();

    /// Return the child block of parameters with given key.
    auto at(const String& key) const -> const Params&;

    /// Return the child parameter with given key.
    auto get(const String& key) const -> const Param&;

    /// Set the child block of parameters with given key to a given block of parameters.
    auto set(const String& key, const Params& node) -> void;

    /// Set the child parameter with given key to a given parameter value.
    auto set(const String& key, const Param& param) -> void;

private:
    Map<String, Any> tree;
};

} // namespace Reaktoro
