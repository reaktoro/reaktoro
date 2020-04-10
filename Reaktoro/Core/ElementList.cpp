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

#include "ElementList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>

namespace Reaktoro {

ElementList::ElementList()
{}

ElementList::ElementList(std::initializer_list<Element> elements)
: m_elements(std::move(elements))
{}

ElementList::ElementList(std::vector<Element> elements)
: m_elements(std::move(elements))
{}

auto ElementList::append(Element element) -> void
{
    m_elements.emplace_back(std::move(element));
}

auto ElementList::data() const -> const std::vector<Element>&
{
    return m_elements;
}

auto ElementList::size() const -> Index
{
    return data().size();
}

auto ElementList::operator[](Index index) const -> const Element&
{
    return data()[index];
}

auto ElementList::indexWithSymbol(String symbol) const -> Index
{
    return indexfn(data(), [&](auto&& e) { return e.symbol() == symbol; });
}

auto ElementList::indexWithName(String name) const -> Index
{
    return indexfn(data(), [&](auto&& e) { return e.name() == name; });
}

auto ElementList::withSymbols(const StringList& symbols) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return contains(symbols, e.symbol()); });
}

auto ElementList::withNames(const StringList& names) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return contains(names, e.name()); });
}

auto ElementList::withTag(String tag) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return contains(e.tags(), tag); });
}

auto ElementList::withoutTag(String tag) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return !contains(e.tags(), tag); });
}

auto ElementList::withTags(const StringList& tags) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return contained(tags, e.tags()); });
}

auto ElementList::withoutTags(const StringList& tags) const -> ElementList
{
    return filter(*this, [&](auto&& e) { return !contained(tags, e.tags()); });
}

} // namespace Reaktoro
