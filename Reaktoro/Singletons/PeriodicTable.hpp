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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

// Forward declarations
class StringList;

/// A type used store a collection of elements.
/// @see Element
class PeriodicTable
{
public:
    /// Construct a copy of a PeriodicTable object [deleted].
    PeriodicTable(const PeriodicTable&) = delete;

    /// Assign a PeriodicTable object to this [deleted].
    auto operator=(const PeriodicTable&) -> PeriodicTable& = delete;

    /// Return the single PeriodicTable object.
    static auto instance() -> PeriodicTable&;

    /// Return the elements in the periodic table.
    static auto elements() -> const Vec<Element>&;

    /// Append a custom element to the periodic table.
    static auto append(Element element) -> void;

    /// Return the number of elements in the periodic table.
    static auto size() -> std::size_t;

    /// Return the element with given symbol.
    static auto elementWithSymbol(String symbol) -> Optional<Element>;

    /// Return the element with given name.
    static auto elementWithName(String name) -> Optional<Element>;

    /// Return the elements with a given tag.
    static auto elementsWithTag(String tag) -> Vec<Element>;

    /// Return the elements with given tags.
    static auto elementsWithTags(const StringList& tags) -> Vec<Element>;

    /// Return begin const iterator of this PeriodicTable instance
    auto begin() const;

    /// Return begin iterator of this PeriodicTable instance
    auto begin();

    /// Return end const iterator of this PeriodicTable instance
    auto end() const;

    /// Return end iterator of this PeriodicTable instance
    auto end();

private:
    /// The elements stored in the periodic table.
    Vec<Element> m_elements;

private:
    /// Construct a default PeriodicTable object [private].
    PeriodicTable();

    /// Destroy this PeriodicTable object [private].
    ~PeriodicTable();
};

} // namespace Reaktoro
