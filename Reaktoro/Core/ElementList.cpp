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

#include "ElementList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>

namespace Reaktoro {

ElementList::ElementList()
{}

ElementList::ElementList(std::initializer_list<Element> elements)
: m_elements(std::move(elements))
{}

ElementList::ElementList(const Vec<Element>& elements)
: m_elements(elements)
{}

auto ElementList::append(const Element& element) -> void
{
    m_elements.push_back(element);
}

auto ElementList::data() const -> const Vec<Element>&
{
    return m_elements;
}

auto ElementList::size() const -> Index
{
    return m_elements.size();
}

auto ElementList::operator[](Index i) const -> const Element&
{
    return m_elements[i];
}

auto ElementList::find(const String& symbol) const -> Index
{
    return findWithSymbol(symbol);
}

auto ElementList::findWithSymbol(const String& symbol) const -> Index
{
    return indexfn(m_elements, RKT_LAMBDA(e, e.symbol() == symbol));
}

auto ElementList::findWithName(const String& name) const -> Index
{
    return indexfn(m_elements, RKT_LAMBDA(e, e.name() == name));
}

auto ElementList::index(const String& symbol) const -> Index
{
    return indexWithSymbol(symbol);
}

auto ElementList::indexWithSymbol(const String& symbol) const -> Index
{
    const auto idx = findWithSymbol(symbol);
    error(idx >= size(), "Could not find any Element object with symbol ", symbol, ".");
    return idx;
}

auto ElementList::indexWithName(const String& name) const -> Index
{
    const auto idx = findWithName(name);
    error(idx >= size(), "Could not find any Element object with name ", name, ".");
    return idx;
}

auto ElementList::get(const String& symbol) const -> const Element&
{
    return getWithSymbol(symbol);
}

auto ElementList::getWithSymbol(const String& symbol) const -> const Element&
{
    return m_elements[indexWithSymbol(symbol)];
}

auto ElementList::getWithName(const String& name) const -> const Element&
{
    return m_elements[indexWithName(name)];
}

auto ElementList::withSymbols(const StringList& symbols) const -> ElementList
{
    return vectorize(symbols, RKT_LAMBDA(symbol, m_elements[indexWithSymbol(symbol)]));
}

auto ElementList::withNames(const StringList& names) const -> ElementList
{
    return vectorize(names, RKT_LAMBDA(name, m_elements[indexWithName(name)]));
}

auto ElementList::withTag(const String& tag) const -> ElementList
{
    return filter(m_elements, RKT_LAMBDA(e, contains(e.tags(), tag)));
}

auto ElementList::withoutTag(const String& tag) const -> ElementList
{
    return filter(m_elements, RKT_LAMBDA(e, !contains(e.tags(), tag)));
}

auto ElementList::withTags(const StringList& tags) const -> ElementList
{
    return filter(m_elements, RKT_LAMBDA(e, contained(tags, e.tags())));
}

auto ElementList::withoutTags(const StringList& tags) const -> ElementList
{
    return filter(m_elements, RKT_LAMBDA(e, !contained(tags, e.tags())));
}

ElementList::operator Vec<Element>&()
{
    return m_elements;
}

ElementList::operator Vec<Element>const&() const
{
    return m_elements;
}

} // namespace Reaktoro
