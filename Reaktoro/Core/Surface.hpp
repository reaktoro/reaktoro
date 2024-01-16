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
#include <Reaktoro/Core/SurfaceAreaModel.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;

/// Used to represent a surface across which chemical reactions take place.
class Surface
{
public:
    /// Construct a default Surface object.
    Surface();

    /// Construct an Surface object with a unique name.
    explicit Surface(String const& name);

    /// Construct an Surface object with a unique name and area model.
    Surface(String const& name, SurfaceAreaModel const& model);

    /// Return a deep copy of this Surface object.
    auto clone() const -> Surface;

    /// Return a duplicate of this Surface object with new name.
    auto withName(String const& name) const -> Surface;

    /// Return a duplicate of this Surface object with new surface area model.
    auto withAreaModel(SurfaceAreaModel const& model) const -> Surface;

    /// Return the unique name of this surface.
    auto name() const -> String const&;

    /// Return the area model of this surface.
    auto areaModel() const -> SurfaceAreaModel const&;

    /// Calculate the area of the surface for given chemical properties of the system (in m2).
    auto area(ChemicalProps const& props) const -> real;

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
