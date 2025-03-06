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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Data.hpp>

namespace Reaktoro {

/// Used to store and retrieve model parameters.
/// @ingroup Core
class Params
{
public:
    /// Construct a default Params object.
    Params();

    /// Construct a Params object with given parameters as Data object.
    Params(Data const& params);

    /// Return a Params object constructed with a given local file.
    /// @warning An exception is thrown if `path` does not point to a valid local parameters file.
    /// @param path The path, including file name, to the local parameters file.
    static auto fromFile(String const& path) -> Params;

    /// Return a Params object constructed with a given embedded file.
    /// @warning An exception is thrown if `path` does not point to a valid embedded parameters file.
    /// @param path The path, including file name, to the embedded parameters file.
    static auto fromEmbeddedFile(String const& path) -> Params;

    /// Return a Params object constructed with given parameters text contents in YAML format.
    /// @param text The text contents of the parameters as a string.
    static auto fromStringYAML(String const& text) -> Params;

    /// Return a Params object constructed with given parameters text contents in JSON format.
    /// @param text The text contents of the parameters as a string.
    static auto fromStringJSON(String const& text) -> Params;

    /// Return parameters with the given file path within the virtual directory of embedded resources (equivalent to `fromEmbeddedFile`).
    static auto embedded(String const& path) -> Params;

    /// Return parameters with the given local file path (equivalent to `fromFile`).
    static auto local(String const& path) -> Params;

    /// Return the underlying Data object of this Params object.
    auto data() const -> Data const&;

    /// Append another Params object into this.
    auto append(Params const& other) -> Params&;

    /// Append a Data object containing parameters into this Params object.
    auto append(Data const& other) -> Params&;

    /// Append another Params object into this.
    auto operator+=(Params const& other) -> Params&;

    /// Return the set of parameters with given name.
    auto operator[](String const& name) const -> Data const&;

private:
    /// The Data object where the model parameters are stored.
    Data m_data;
};

} // namespace Reaktoro
