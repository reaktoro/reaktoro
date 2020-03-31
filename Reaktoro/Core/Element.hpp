// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <memory>
#include <string>
#include <vector>

namespace Reaktoro {

/// A type used to define a element and its attributes.
class Element
{
public:
    /// A type used to represent the arguments to construct an Element object.
    struct Args
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
    };

    /// Construct a default Element object.
    Element();

    /// Construct an Element object with given symbol.
    explicit Element(std::string symbol);

    /// Construct an Element object with given data.
    /// @param args The arguments to construct the element.
    explicit Element(const Args& args);

    /// Return a duplicate of this Element object with replaced symbol attribute.
    auto withSymbol(std::string symbol) const -> Element;

    /// Return a duplicate of this Element object with replaced name attribute.
    auto withName(std::string name) const -> Element;

    /// Return a duplicate of this Element object with replaced atomic number attribute.
    auto withAtomicNumber(std::size_t atomic_number) const -> Element;

    /// Return a duplicate of this Element object with replaced atomic weight attribute.
    /// @note This method is equivalent to withMolarMass
    auto withAtomicWeight(double atomic_weight) const -> Element;

    /// Return a duplicate of this Element object with replaced molar mass attribute.
    /// @note This method is equivalent to withAtomicWeight
    auto withMolarMass(double molar_mass) const -> Element;

    /// Return a duplicate of this Element object with replaced electronegativity attribute.
    auto withElectronegativity(double electronegativity) const -> Element;

    /// Return a duplicate of this Element object with replaced tags attribute.
    auto withTags(std::vector<std::string> tags) const -> Element;

    /// Return the symbol of the element (e.g., "H", "O", "C", "Na").
    auto symbol() const -> std::string;

    /// Return the name of the element (e.g., "Hydrogen", "Oxygen").
    auto name() const -> std::string;

    /// Return the atomic number of the element.
    auto atomicNumber() const -> std::size_t;

    /// Return the atomic weight of the element (in unit of kg/mol).
    /// @note This method is equivalent to molarMass
    auto atomicWeight() const -> double;

    /// Return the molar mass of the element (in unit of kg/mol).
    /// @note This method is equivalent to atomicWeight
    auto molarMass() const -> double;

    /// Return the electronegativity of the element.
    auto electronegativity() const -> double;

    /// Return the tags of the element.
    auto tags() const -> const std::vector<std::string>&;

    /// Return a deep copy of this Element object.
    auto clone() const -> Element;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
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
