// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

/// Used to represent the interface surface between two phases in a chemical system.
class Surface
{
public:
    /// Construct a default Surface object.
    Surface();

    /// Construct an Surface object with a unique name.
    explicit Surface(String const& name);

    /// Construct an Surface object with a unique name and names of phases between which the surface exists.
    Surface(String const& name, String const& phase1, String const& phase2);

    /// Return a deep copy of this Surface object.
    auto clone() const -> Surface;

    /// Return a duplicate of this Surface object with replaced name.
    auto withName(String const& name) const -> Surface;

    /// Return a duplicate of this Surface object with replaced names of the phases between which this surface exists.
    auto withPhases(String const& phase1, String const& phase2) const -> Surface;

    /// Return the unique name of this surface.
    auto name() const -> String;

    /// Return the names of the phases between which this surface exists.
    auto phases() const -> Pair<String, String> const&;

    /// Return true if this surface is equivalent to another with the same interface phases.
    auto equivalent(Surface const& another) const -> bool;

    /// Return true if this surface is equivalent to another with the same interface phases using phase names.
    auto equivalent(String const& phase1, String const& phase2) const -> bool;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Compare two Surface objects for less than.
auto operator<(Surface const& lhs, Surface const& rhs) -> bool;

/// Compare two Surface objects for equality.
auto operator==(Surface const& lhs, Surface const& rhs) -> bool;

} // namespace Reaktoro

// Custom specialization of std::hash for Reaktoro::Surface
namespace std {

template<>
struct hash<Reaktoro::Surface>
{
    std::size_t operator()(Reaktoro::Surface const& s) const noexcept
    {
        return std::hash<std::string>{}(s.name());
    }
};

} // namespace std
