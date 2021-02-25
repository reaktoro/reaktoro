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

#include "Element.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Singletons/PeriodicTable.hpp>

namespace Reaktoro {
namespace detail {

/// Return an Element object from PeriodicTable with given symbol.
auto getElementFromPeriodicTable(String symbol) -> Element
{
    const auto element = PeriodicTable::elementWithSymbol(symbol);
    error(!element.has_value(), "Cannot proceed with Element(symbol) constructor. "
        "PeriodicTable contains no element with symbol ", symbol, ". "
        "Use method PeriodicTable::append (in C++) or PeriodicTable.append (in Python) "
        "to add a new Element with this symbol. Or construct the Element object without "
        "using this special constructor.");
    return element.value();
}

} // namespace detail

struct Element::Impl
{
    /// The symbol of the element (e.g., "H", "O", "C", "Na").
    String symbol;

    /// The name of the element (e.g., "Hydrogen", "Oxygen").
    String name;

    /// The atomic number of the element.
    Index atomic_number = {};

    /// The atomic weight (or molar mass) of the element (in kg/mol).
    double atomic_weight = {};

    /// The electronegativity of the element.
    double electronegativity = {};

    /// The tags of the element.
    Strings tags;

    /// Construct a default Element::Impl object.
    Impl()
    {}

    /// Construct an Element::Impl object with given attributes.
    Impl(const Args& args)
    : symbol(args.symbol),
      name(args.name),
      atomic_number(args.atomic_number),
      atomic_weight(args.atomic_weight),
      electronegativity(args.electronegativity),
      tags(args.tags)
    {}
};

Element::Element()
: pimpl(new Impl())
{}

Element::Element(String symbol)
 : Element(detail::getElementFromPeriodicTable(symbol))
{}

Element::Element(const Args& attributes)
 : pimpl(new Impl(attributes))
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

auto Element::withName(String name) const -> Element
{
    Element copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Element::withAtomicNumber(Index atomic_number) const -> Element
{
    Element copy = clone();
    copy.pimpl->atomic_number = atomic_number;
    return copy;
}

auto Element::withAtomicWeight(double atomic_weight) const -> Element
{
    Element copy = clone();
    copy.pimpl->atomic_weight = atomic_weight;
    return copy;
}

auto Element::withMolarMass(double molar_mass) const -> Element
{
    return withAtomicWeight(molar_mass);
}

auto Element::withElectronegativity(double electronegativity) const -> Element
{
    Element copy = clone();
    copy.pimpl->electronegativity = electronegativity;
    return copy;
}

auto Element::withTags(Strings tags) const -> Element
{
    Element copy = clone();
    copy.pimpl->tags = tags;
    return copy;
}

auto Element::symbol() const -> String
{
    return pimpl->symbol;
}

auto Element::name() const -> String
{
    return pimpl->name;
}

auto Element::atomicNumber() const -> Index
{
    return pimpl->atomic_number;
}

auto Element::atomicWeight() const -> double
{
    return pimpl->atomic_weight;
}

auto Element::molarMass() const -> double
{
    return atomicWeight();
}

auto Element::electronegativity() const -> double
{
    return pimpl->electronegativity;
}

auto Element::tags() const -> const Strings&
{
    return pimpl->tags;
}

auto operator<(const Element& lhs, const Element& rhs) -> bool
{
    return lhs.atomicNumber() < rhs.atomicNumber();
}

auto operator==(const Element& lhs, const Element& rhs) -> bool
{
    return lhs.symbol()            == rhs.symbol()            &&
           lhs.name()              == rhs.name()              &&
           lhs.atomicNumber()      == rhs.atomicNumber()      &&
           lhs.atomicWeight()      == rhs.atomicWeight()      &&
           lhs.electronegativity() == rhs.electronegativity() &&
           lhs.tags()              == rhs.tags()
           ;
}

} // namespace Reaktoro
