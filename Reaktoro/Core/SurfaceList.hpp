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
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/Surface.hpp>

namespace Reaktoro {

/// A type used as a collection of surfaces.
class SurfaceList
{
public:
    /// Construct a default SurfaceList object.
    SurfaceList();

    /// Construct an SurfaceList object with given surfaces.
    SurfaceList(std::initializer_list<Surface> surfaces);

    /// Construct an SurfaceList object with given surfaces.
    SurfaceList(Vec<Surface> const& surfaces);

    /// Append a new surface to the list of surfaces.
    auto append(Surface const& surface) -> void;

    /// Return the internal collection of Surface objects.
    auto data() const -> Vec<Surface> const&;

    /// Return true if there are no surfaces in the collection.
    auto empty() const -> bool;

    /// Return the number of surfaces in the collection.
    auto size() const -> Index;

    /// Return the Surface object with given index.
    auto operator[](Index i) const -> Surface const&;

    /// Return the Surface object with given index.
    auto operator[](Index i) -> Surface&;

    /// Return the index of the first surface with given name or the number of surfaces if not found.
    auto find(String const& name) const -> Index;

    /// Return the index of the first surface with given unique name or the number of surfaces if not found.
    auto findWithName(String const& name) const -> Index;

    /// Return the index of the first surface with given name or throw a runtime error if not found.
    auto index(String const& name) const -> Index;

    /// Return the index of the first surface with given unique name or throw a runtime error if not found.
    auto indexWithName(String const& name) const -> Index;

    /// Return the surface with given name.
    auto get(String const& name) const -> Surface const&;

    /// Return the surface with given name or throw a runtime error if not found.
    auto getWithName(String const& name) const -> Surface const&;

    /// Return all surfaces with given names.
    auto withNames(StringList const& names) const -> SurfaceList;

    /// Convert this SurfaceList object into its Vec<Surface>.
    operator Vec<Surface>&();

    /// Convert this SurfaceList object into its Vec<Surface>.
    operator Vec<Surface>const&() const;

private:
    /// The surfaces stored in the list.
    Vec<Surface> m_surfaces;

public:
    /// Construct an SurfaceList object with given begin and end iterators.
    template<typename InputIterator>
    SurfaceList(InputIterator begin, InputIterator end) : m_surfaces(begin, end) {}

    /// Return begin const iterator of this SurfaceList instance (for STL compatibility reasons).
    auto begin() const { return m_surfaces.begin(); }

    /// Return begin iterator of this SurfaceList instance (for STL compatibility reasons).
    auto begin() { return m_surfaces.begin(); }

    /// Return end const iterator of this SurfaceList instance (for STL compatibility reasons).
    auto end() const { return m_surfaces.end(); }

    /// Return end iterator of this SurfaceList instance (for STL compatibility reasons).
    auto end() { return m_surfaces.end(); }

    /// Append a new Surface at the back of the container (for STL compatibility reasons).
    auto push_back(Surface const& surfaces) -> void { append(surfaces); }

    /// Insert a container of Surface objects into this SurfaceList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_surfaces.insert(pos, begin, end); }

    /// The type of the value stored in a SurfaceList (for STL compatibility reasons).
    using value_type = Surface;
};

/// Return the concatenation of two SurfaceList objects.
auto operator+(SurfaceList const& a, SurfaceList const& b) -> SurfaceList;

} // namespace Reaktoro
