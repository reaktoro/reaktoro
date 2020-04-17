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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

// Forward declaration of ElementListBase
template<typename Data>
class ElementListBase;

/// The specialized container to deal with a collection of Element objects.
using ElementList = ElementListBase<Vec<Element>>;

/// The specialized container to deal with a const reference view of a collection of Element objects.
using ElementListConstRef = ElementListBase<const Vec<Element>&>;

/// A type used as a collection of elements.
template<typename Data>
class ElementListBase
{
public:
    /// Construct a default ElementListBase object.
    ElementListBase() {}

    /// Construct an ElementListBase object with given elements.
    ElementListBase(std::initializer_list<Element> elements) : m_elements(std::move(elements)) {}

    /// Construct an ElementListBase object with given elements.
    ElementListBase(const Vec<Element>& elements) : m_elements(elements) {}

    /// Construct an ElementListBase object with given another one.
    template<typename OtherData>
    ElementListBase(const ElementListBase<OtherData>& other) : m_elements(other.m_elements) {}

    /// Append a new element to the list of elements.
    auto append(const Element& element)
    {
        m_elements.push_back(element);
    }

    /// Return the internal collection of Element objects.
    auto data() const -> const Vec<Element>&
    {
        return m_elements;
    }

    /// Return the number of elements in the collection.
    auto size() const -> Index
    {
        return m_elements.size();
    }

    /// Return the Element object with given index.
    auto operator[](Index i) const -> const Element&
    {
        return m_elements[i];
    }

    /// Return the index of the first element with given symbol or number of elements if not found.
    auto indexWithSymbol(const String& symbol) const -> Index
    {
        return indexfn(m_elements, RKT_LAMBDA(e, e.symbol() == symbol));
    }

    /// Return the index of the first element with given unique name or number of elements if not found.
    auto indexWithName(const String& name) const -> Index
    {
        return indexfn(m_elements, RKT_LAMBDA(e, e.name() == name));
    }

    /// Return all elements with given symbols.
    auto withSymbols(const StringList& symbols) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, contains(symbols, e.symbol())));
    }

    /// Return all elements with given names.
    auto withNames(const StringList& names) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, contains(names, e.name())));
    }

    /// Return all elements with a given tag.
    auto withTag(const String& tag) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, contains(e.tags(), tag)));
    }

    /// Return all elements without a given tag.
    auto withoutTag(const String& tag) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, !contains(e.tags(), tag)));
    }

    /// Return all elements with given tags.
    auto withTags(const StringList& tags) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, contained(tags, e.tags())));
    }

    /// Return all elements without given tags.
    auto withoutTags(const StringList& tags) const -> ElementList
    {
        return filter(m_elements, RKT_LAMBDA(e, !contained(tags, e.tags())));
    }

private:
    /// The elements stored in the list.
    Vec<Element> m_elements;

public:
    /// Construct an ElementListBase object with given begin and end iterators.
    template<typename InputIterator>
    ElementListBase(InputIterator begin, InputIterator end) : m_elements(begin, end) {}

    /// Return begin const iterator of this ElementListBase instance (for STL compatibility reasons).
    inline auto begin() const { return data().begin(); }

    /// Return begin iterator of this ElementListBase instance (for STL compatibility reasons).
    inline auto begin() { return data().begin(); }

    /// Return end const iterator of this ElementListBase instance (for STL compatibility reasons).
    inline auto end() const { return data().end(); }

    /// Return end iterator of this ElementListBase instance (for STL compatibility reasons).
    inline auto end() { return data().end(); }

    /// Append a new Element at the back of the container (for STL compatibility reasons).
    inline auto push_back(const Element& elements) -> void { append(elements); }

    /// The type of the value stored in a ElementListBase (for STL compatibility reasons).
    using value_type = Element;
};

} // namespace Reaktoro
