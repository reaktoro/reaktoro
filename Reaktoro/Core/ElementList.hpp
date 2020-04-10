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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

/// A type used as a collection of elements.
class ElementList
{
public:
    /// Construct a default ElementList object.
    ElementList();

    /// Construct an ElementList object with given elements.
    ElementList(std::initializer_list<Element> elements);

    /// Construct an ElementList object with given elements.
    explicit ElementList(std::vector<Element> elements);

    /// Append a new element to the list of elements.
    auto append(Element element) -> void;

    /// Return the internal collection of Element objects.
    auto data() const -> const std::vector<Element>&;

    /// Return the number of elements in the collection.
    auto size() const -> Index;

    /// Return the Element object with given index.
    auto operator[](Index index) const -> const Element&;

    /// Return the index of the first element with equivalent substance formula.
    /// If there is no elements with given substance formula, return number of elements.
    auto indexWithSymbol(String formula) const -> Index;

    /// Return the index of the first element with given unique name.
    /// If there is no elements with given name, return number of elements.
    auto indexWithName(String name) const -> Index;

    /// Return all elements with given symbols.
    auto withSymbols(const StringList& symbols) const -> ElementList;

    /// Return all elements with given names.
    auto withNames(const StringList& names) const -> ElementList;

    /// Return all elements with a given tag.
    auto withTag(String tag) const -> ElementList;

    /// Return all elements without a given tag.
    auto withoutTag(String tag) const -> ElementList;

    /// Return all elements with given tags.
    auto withTags(const StringList& tags) const -> ElementList;

    /// Return all elements without given tags.
    auto withoutTags(const StringList& tags) const -> ElementList;

private:
    /// The elements stored in the list.
    std::vector<Element> m_elements;

public:
    /// Construct an ElementList object with given begin and end iterators.
    template<typename InputIterator>
    ElementList(InputIterator begin, InputIterator end) : m_elements(begin, end) {}

    /// Return begin const iterator of this ElementList instance (for STL compatibility reasons).
    inline auto begin() const { return data().begin(); }

    /// Return begin iterator of this ElementList instance (for STL compatibility reasons).
    inline auto begin() { return data().begin(); }

    /// Return end const iterator of this ElementList instance (for STL compatibility reasons).
    inline auto end() const { return data().end(); }

    /// Return end iterator of this ElementList instance (for STL compatibility reasons).
    inline auto end() { return data().end(); }

    /// Append a new Element at the back of the container (for STL compatibility reasons).
    inline auto push_back(const Element& elements) -> void { append(elements); }

    /// The type of the value stored in a ElementList (for STL compatibility reasons).
    using value_type = Element;
};

} // namespace Reaktoro
