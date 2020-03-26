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
#include <any>
#include <map>
#include <memory>
#include <string>

namespace Reaktoro {

// Forward declarations
class Element;

/// A type used to describe a species and its attributes.
/// The Species class is used to represent a species. It is an important
/// class in the library, since it defines fundamental attributes of a general
/// species such as its elemental formula, electrical charge and molar mass.
/// @see Phase
/// @ingroup Core
class Species
{
public:
    /// Construct a default Species instance.
    Species();

    /// Return a copy of this Species object with given name.
    auto withName(std::string name) -> Species;

    /// Return a copy of this Species object with given formula.
    auto withFormula(std::string formula) -> Species;

    /// Return a copy of this Species object with given elements.
    auto withElements(const std::map<Element, double>& elements) -> Species;

    /// Return a copy of this Species object with given type.
    /// The following are examples of species types:
    /// `Aqueous`, `Gaseous`, `Liquid`, `Mineral`.
    auto withType(std::string type) -> Species;

    /// Return a copy of this Species object with given data.
    auto withData(const std::any& data) -> Species;

    /// Return the number of elements of the species.
    auto numElements() const -> unsigned;

    /// Return the name of the species.
    auto name() const -> std::string;

    /// Return the formula of the species.
    auto formula() const -> std::string;

    /// Return the elements that compose the species and their coefficients.
    auto elements() const -> const std::map<Element, double>&;

    /// Return the type of the species.
    auto type() const -> std::string;

    /// Return the molar mass of the species (in units of kg/mol).
    auto molarMass() const -> double;

    /// Return the electrical charge of the species.
    auto charge() const -> double;

    /// Return the stoichiometry of an element in the species.
    auto elementCoefficient(std::string element) const -> double;

    /// Return the specific data the species may have such as thermodynamic data.
    auto data() const -> const std::any&;

    /// Return a cloned copy of this Species object.
    auto clone() const -> Species;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Species instances for less than
auto operator<(const Species& lhs, const Species& rhs) -> bool;

/// Compare two Species instances for equality
auto operator==(const Species& lhs, const Species& rhs) -> bool;

} // namespace Reaktoro
