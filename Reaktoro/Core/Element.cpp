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

#include "Element.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>

namespace Reaktoro {

struct Element::Impl
{
    /// The symbol of the element (e.g., "H", "O", "C", "Na").
    std::string symbol;

    /// The name of the element (e.g., "Hydrogen", "Oxygen").
    std::string name;

    /// The atomic number of the element.
    std::size_t atomic_number = {};

    /// The atomic weight (or molar mass) of the element (in unit of kg/mol).
    double atomic_weight = {};

    /// The electronegativity of the element.
    double electronegativity = {};

    /// The tags of the element.
    std::vector<std::string> tags;

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

Element::Element(std::string symbol)
 : pimpl(new Impl({symbol}))
{}

Element::Element(const Args& attributes)
 : pimpl(new Impl(attributes))
{}

auto Element::withSymbol(std::string symbol) const -> Element
{
    Element copy = clone();
    copy.pimpl->symbol = symbol;
    return copy;
}

auto Element::withName(std::string name) const -> Element
{
    Element copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Element::withAtomicNumber(std::size_t atomic_number) const -> Element
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

auto Element::withTags(std::vector<std::string> tags) const -> Element
{
    Element copy = clone();
    copy.pimpl->tags = tags;
    return copy;
}

auto Element::symbol() const -> std::string
{
    return pimpl->symbol;
}

auto Element::name() const -> std::string
{
    return pimpl->name;
}

auto Element::atomicNumber() const -> std::size_t
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

auto Element::tags() const -> const std::vector<std::string>&
{
    return pimpl->tags;
}

auto Element::clone() const -> Element
{
    Element element;
    *element.pimpl = *pimpl;
    return element;
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
