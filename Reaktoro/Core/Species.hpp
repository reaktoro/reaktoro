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
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Reaktoro {

// Forward declarations (class)
class ChemicalFormula;
class Element;

// Forward declarations (enum)
enum class AggregateState;

/// A type used to represent a chemical species and its attributes.
class Species
{
public:
    /// Construct a default Species object.
    Species();

    /// Construct a Species object with given chemical formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    Species(const ChemicalFormula& formula);

    /// Return a duplicate of this Species object with replaced name attribute.
    auto withName(std::string name) -> Species;

    /// Return a duplicate of this Species object with replaced formula attribute.
    auto withFormula(const ChemicalFormula& formula) -> Species;

    /// Return a duplicate of this Species object with replaced aggregate state.
    auto withAggregateState(AggregateState option) -> Species;

    /// Return a duplicate of this Species object with replaced tags attribute.
    auto withTags(std::vector<std::string> tags) -> Species;

    /// Return a duplicate of this Species object with attached data.
    auto withAttachedData(std::string id, std::any data) -> Species;

    /// Return the name of the species if provided, otherwise, its formula.
    auto name() const -> std::string;

    /// Return the chemical formula of the species.
    auto formula() const -> const ChemicalFormula&;

    /// Return the electric charge of the species.
    auto charge() const -> double;

    /// Return the molar mass of the species (in unit of kg/mol).
    auto molarMass() const -> double;

    /// Return the aggregate state of the species.
    auto aggregateState() const -> AggregateState;

    /// Return the elements that compose the species and their coefficients.
    auto elements() const -> const std::vector<std::pair<Element, double>>&;

    /// Return the coefficient of an element in the species.
    auto elementCoefficient(const std::string& symbol) const -> double;

    /// Return the tags of the species (e.g., `organic`, `mineral`).
    auto tags() const -> const std::vector<std::string>&;

    /// Return the attached data with given id.
    auto attachedData(std::string id) const -> std::optional<std::any>;

    /// Return all attached data and their corresponding id's.
    auto attachedData() const -> const std::unordered_map<std::string, std::any>&;

    /// Return a deep copy of this Species object.
    auto clone() const -> Species;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Species objects for less than
auto operator<(const Species& lhs, const Species& rhs) -> bool;

/// Compare two Species objects for equality
auto operator==(const Species& lhs, const Species& rhs) -> bool;

} // namespace Reaktoro
