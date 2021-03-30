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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

// Forward declarations
class yaml;

/// The class used to parse databases and construct a Database object.
/// @ingroup Core
class DatabaseParser
{
public:
    /// Construct a default DatabaseParser object.
    DatabaseParser();

    /// Construct a copy of a DatabaseParser object.
    DatabaseParser(const DatabaseParser& other);

    /// Construct a DatabaseParser object with given YAML node.
    explicit DatabaseParser(const yaml& node);

    /// Destroy this DatabaseParser object.
    ~DatabaseParser();

    /// Assign another DatabaseParser object to this.
    auto operator=(DatabaseParser other) -> DatabaseParser&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
