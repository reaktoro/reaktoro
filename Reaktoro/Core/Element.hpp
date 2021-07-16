// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// A type used to define a element and its attributes.
class Element
{
public:
    /// The attributes of an Element object.
    struct Attribs
    {
        /// The symbol of the element (e.g., "H", "O", "C", "Na").
        /// @note This is a required attribute. The element symbol must be unique among all other elements.
        String symbol;

        /// The molar mass of the element (in kg/mol).
        /// @note This is a required attribute. Zero is allowed, but this
        /// prevents conversion from amount to mass values for chemical species
        /// containing elements with zero molar mass.
        double molar_mass = 0.0;

        /// The name of the element (e.g., "Hydrogen", "Oxygen").
        /// @note This is an optional attribute. If not provided, assumed equal to `symbol`.
        String name;

        /// The tags of the element.
        /// @note This is an optional attribute.
        Strings tags;
    };

    /// Construct a default Element object.
    Element();

    /// Construct an Element object by looking up to the periodic table with given symbol.
    explicit Element(String symbol);

    /// Construct an Element object with given attributes.
    explicit Element(const Attribs& attribs);

    /// Return a deep copy of this Element object.
    auto clone() const -> Element;

    /// Return a duplicate of this Element object with replaced symbol attribute.
    auto withSymbol(String symbol) const -> Element;

    /// Return a duplicate of this Element object with replaced molar mass attribute (in kg/mol).
    auto withMolarMass(double value) const -> Element;

    /// Return a duplicate of this Element object with replaced name attribute.
    auto withName(String name) const -> Element;

    /// Return a duplicate of this Element object with replaced tags attribute.
    auto withTags(const StringList& tags) const -> Element;

    /// Return the symbol of the element (e.g., "H", "O", "C", "Na").
    auto symbol() const -> String;

    /// Return the molar mass of the element (in kg/mol).
    auto molarMass() const -> double;

    /// Return the name of the element (e.g., "Hydrogen", "Oxygen").
    auto name() const -> String;

    /// Return the tags of the element.
    auto tags() const -> const Strings&;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Compare two Element objects for less than.
auto operator<(const Element& lhs, const Element& rhs) -> bool;

/// Compare two Element objects for equality.
auto operator==(const Element& lhs, const Element& rhs) -> bool;

} // namespace Reaktoro

// Custom specialization of std::hash for Reaktoro::Element
namespace std {

template<>
struct hash<Reaktoro::Element>
{
    std::size_t operator()(const Reaktoro::Element& e) const noexcept
    {
        return std::hash<std::string>{}(e.symbol());
    }
};

} // namespace std
