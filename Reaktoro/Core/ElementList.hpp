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
    ElementList(const Vec<Element>& elements);

    /// Append a new element to the list of elements.
    auto append(const Element& element) -> void;

    /// Return the internal collection of Element objects.
    auto data() const -> const Vec<Element>&;

    /// Return true if there are no elements in the collection.
    auto empty() const -> bool;

    /// Return the number of elements in the collection.
    auto size() const -> Index;

    /// Return the Element object with given index.
    auto operator[](Index i) const -> const Element&;

    /// Return the index of the first element with given symbol or the number of elements if not found.
    auto find(const String& symbol) const -> Index;

    /// Return the index of the first element with given symbol or the number of elements if not found.
    auto findWithSymbol(const String& symbol) const -> Index;

    /// Return the index of the first element with given unique name or the number of elements if not found.
    auto findWithName(const String& name) const -> Index;

    /// Return the index of the first element with given symbol or throw a runtime error if not found.
    auto index(const String& symbol) const -> Index;

    /// Return the index of the first element with given symbol or throw a runtime error if not found.
    auto indexWithSymbol(const String& symbol) const -> Index;

    /// Return the index of the first element with given unique name or throw a runtime error if not found.
    auto indexWithName(const String& name) const -> Index;

    /// Return the element with given symbol.
    auto get(const String& symbol) const -> const Element&;

    /// Return the element with given symbol.
    auto getWithSymbol(const String& symbol) const -> const Element&;

    /// Return the element with given name.
    auto getWithName(const String& name) const -> const Element&;

    /// Return all elements with given symbols.
    auto withSymbols(const StringList& symbols) const -> ElementList;

    /// Return all elements with given names.
    auto withNames(const StringList& names) const -> ElementList;

    /// Return all elements with a given tag.
    auto withTag(const String& tag) const -> ElementList;

    /// Return all elements without a given tag.
    auto withoutTag(const String& tag) const -> ElementList;

    /// Return all elements with given tags.
    auto withTags(const StringList& tags) const -> ElementList;

    /// Return all elements without given tags.
    auto withoutTags(const StringList& tags) const -> ElementList;

    /// Convert this ElementList object into its Vec<Element>.
    operator Vec<Element>&();

    /// Convert this ElementList object into its Vec<Element>.
    operator Vec<Element>const&() const;

private:
    /// The elements stored in the list.
    Vec<Element> m_elements;

public:
    /// Construct an ElementList object with given begin and end iterators.
    template<typename InputIterator>
    ElementList(InputIterator begin, InputIterator end) : m_elements(begin, end) {}

    /// Return begin const iterator of this ElementList instance (for STL compatibility reasons).
    auto begin() const { return m_elements.begin(); }

    /// Return begin iterator of this ElementList instance (for STL compatibility reasons).
    auto begin() { return m_elements.begin(); }

    /// Return end const iterator of this ElementList instance (for STL compatibility reasons).
    auto end() const { return m_elements.end(); }

    /// Return end iterator of this ElementList instance (for STL compatibility reasons).
    auto end() { return m_elements.end(); }

    /// Append a new Element at the back of the container (for STL compatibility reasons).
    auto push_back(const Element& elements) -> void { append(elements); }

    /// Insert a container of Element objects into this ElementList instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_elements.insert(pos, begin, end); }

    /// The type of the value stored in a ElementList (for STL compatibility reasons).
    using value_type = Element;
};

/// Return the concatenation of two ElementList objects.
auto operator+(const ElementList& a, const ElementList& b) -> ElementList;

} // namespace Reaktoro
