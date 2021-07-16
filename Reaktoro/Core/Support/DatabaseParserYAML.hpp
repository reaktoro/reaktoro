// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

// Forward declarations
class Database;
class yaml;

/// The auxiliary class used to handle the parsing of YAML files to construct a Database object.
class DatabaseParserYAML
{
public:
    /// Construct a default DatabaseParserYAML object.
    DatabaseParserYAML();

    /// Construct a copy of a DatabaseParserYAML object.
    DatabaseParserYAML(const DatabaseParserYAML& other);

    /// Construct a DatabaseParserYAML object with given YAML node.
    explicit DatabaseParserYAML(const yaml& node);

    /// Destroy this DatabaseParserYAML object.
    ~DatabaseParserYAML();

    /// Assign another DatabaseParserYAML object to this.
    auto operator=(DatabaseParserYAML other) -> DatabaseParserYAML&;

    /// Return the parsed Element objects in the database file.
    auto elements() const -> const ElementList&;

    /// Return the parsed Species objects in the database file.
    auto species() const -> const SpeciesList&;

    /// Return the parsed Element objects in the database file.
    operator Database() const;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
