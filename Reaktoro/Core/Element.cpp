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

#include "Element.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Singletons/Elements.hpp>

namespace Reaktoro {
namespace detail {

/// Return an Element object from Elements with given symbol.
auto getDefaultElement(String symbol) -> Element
{
    const auto element = Elements::withSymbol(symbol);
    errorif(!element.has_value(), "Cannot proceed with Element(symbol) constructor. "
        "Could not find an element in Elements with symbol ", symbol, ". "
        "Use method Elements::append (in C++) or Elements.append (in Python) "
        "to add a new Element with this symbol. Or construct the Element object without "
        "using this special constructor.");
    return element.value();
}

} // namespace detail

struct Element::Impl
{
    /// The symbol of the element (e.g., "H", "O", "C", "Na").
    String symbol;

    /// The molar mass of the element (in kg/mol).
    double molar_mass = {};

    /// The name of the element (e.g., "Hydrogen", "Oxygen").
    String name;

    /// The tags of the element.
    Strings tags;

    /// Construct a default Element::Impl object.
    Impl()
    {}

    /// Construct an Element::Impl object with given attributes.
    Impl(const Attribs& attribs)
    {
        errorif(attribs.symbol.empty(), "Could not construct Element object with constructor Element(Element::Attribs). "
            "Element::Attribs::symbol cannot be empty.");
        errorif(attribs.molar_mass < 0.0, "Could not construct Element object with constructor Element(Element::Attribs). "
            "Element::Attribs::molar_mass cannot be negative.");
        symbol = attribs.symbol;
        molar_mass = attribs.molar_mass;
        name = attribs.name.empty() ? symbol : attribs.name;
        tags = attribs.tags;
    }
};

Element::Element()
: pimpl(new Impl())
{}

Element::Element(String symbol)
 : Element(detail::getDefaultElement(symbol))
{}

Element::Element(const Attribs& attribs)
 : pimpl(new Impl(attribs))
{}

auto Element::clone() const -> Element
{
    Element element;
    *element.pimpl = *pimpl;
    return element;
}

auto Element::withSymbol(String symbol) const -> Element
{
    Element copy = clone();
    copy.pimpl->symbol = symbol;
    return copy;
}

auto Element::withMolarMass(double molar_mass) const -> Element
{
    errorif(molar_mass < 0.0, "Cannot set the molar mass of an Element to a negative value");
    Element copy = clone();
    copy.pimpl->molar_mass = molar_mass;
    return copy;
}

auto Element::withName(String name) const -> Element
{
    Element copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Element::withTags(const StringList& tags) const -> Element
{
    Element copy = clone();
    copy.pimpl->tags = tags;
    return copy;
}

auto Element::symbol() const -> String
{
    return pimpl->symbol;
}

auto Element::molarMass() const -> double
{
    return pimpl->molar_mass;
}

auto Element::name() const -> String
{
    return pimpl->name;
}

auto Element::tags() const -> const Strings&
{
    return pimpl->tags;
}

auto operator<(const Element& lhs, const Element& rhs) -> bool
{
    return lhs.symbol() < rhs.symbol();
}

auto operator==(const Element& lhs, const Element& rhs) -> bool
{
    return lhs.symbol()            == rhs.symbol()            &&
           lhs.molarMass()         == rhs.molarMass()         &&
           lhs.name()              == rhs.name()              &&
           lhs.tags()              == rhs.tags()
           ;
}

} // namespace Reaktoro
