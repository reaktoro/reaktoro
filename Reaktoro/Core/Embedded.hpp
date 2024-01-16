// Reaktoro is a unified framework for modeling chemically reactive phases.
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// Used to retrieve embedded resources (e.g., database files, parameter files) in Reaktoro.
class Embedded
{
public:
    /// Return the contents of the embedded document with given path (as a string).
    static auto get(String const& path) -> String;

    /// Return the contents of the embedded document with given path (as a string).
    static auto getAsString(String const& path) -> String;

    /// Return the contents of the embedded document with given path (as a string view).
    static auto getAsStringView(String const& path) -> Pair<Chars, Chars>;

    /// Deleted default constructor.
    Embedded() = delete;
};

} // namespace Reaktoro
