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

namespace Reaktoro {

// Forward declarations
class Database;
class yaml;

/// The auxiliary class used to construct a Database object from an YAML document (decoder).
class DatabaseDecoderYAML
{
public:
    /// Construct a default DatabaseDecoderYAML object.
    DatabaseDecoderYAML();

    /// Construct a copy of a DatabaseDecoderYAML object.
    DatabaseDecoderYAML(const DatabaseDecoderYAML& other);

    /// Construct a DatabaseDecoderYAML object with given YAML node.
    explicit DatabaseDecoderYAML(const yaml& node);

    /// Destroy this DatabaseDecoderYAML object.
    ~DatabaseDecoderYAML();

    /// Assign another DatabaseDecoderYAML object to this.
    auto operator=(DatabaseDecoderYAML other) -> DatabaseDecoderYAML&;

    /// Convert this DatabaseDecoderYAML object into a Database object.
    operator Database() const;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
