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
#include <string>
#include <unordered_map>
#include <vector>

namespace Reaktoro {

// Forward declarations (class)
class Element;
class SubstanceCriticalProps;

// Forward declarations (enum)
enum class AggregateState;

/// A type used to represent a chemical species and its attributes.
class Species
{
public:
    /// The type for the container of Element objects composing the species and their corresponding coefficients.
    using Elements = std::unordered_map<Element, double>;

    /// The type for the container of chemical element symbols composing the species and their corresponding coefficients.
    using ElementSymbols = std::unordered_map<std::string, double>;

    /// Construct a default Species object.
    Species();

    /// Construct a Species object with given chemical formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    Species(std::string formula);

    /// Return a duplicate of this Species object with new symbol that uniquely identifies this species.
    auto withSymbol(std::string symbol) -> Species;

    /// Return a duplicate of this Species object with new substance name attribute (case insensitive).
    auto withName(std::string name) -> Species;

    /// Return a duplicate of this Species object with new formula attribute.
    auto withFormula(std::string formula) -> Species;

    /// Return a duplicate of this Species object with new Element objects and respective coefficients.
    auto withElements(Elements elements) -> Species;

    /// Return a duplicate of this Species object with new Element objects and respective coefficients.
    /// @note The Element objects in the Species instance will be collected from PeriodicTable.
    auto withElementSymbols(ElementSymbols symbols) -> Species;

    /// Return a duplicate of this Species object with new electric charge attribute.
    auto withCharge(double charge) -> Species;

    /// Return a duplicate of this Species object with new aggregate state.
    auto withAggregateState(AggregateState option) -> Species;

    /// Return a duplicate of this Species object with new tags attribute.
    auto withTags(std::vector<std::string> tags) -> Species;

    /// Return a duplicate of this Species object with new critical properties.
    auto withCriticalProps(const SubstanceCriticalProps& props) -> Species;

    /// Return a duplicate of this Species object with new attached data whose type is known at runtime only.
    auto withAttachedData(std::any data) -> Species;

    /// Return the symbol that uniquely identifies this species if provided, otherwise, its formula.
    auto symbol() const -> std::string;

    /// Return the name of the underlying substance of the species if provided, otherwise, its formula.
    auto name() const -> std::string;

    /// Return the chemical formula of the species.
    auto formula() const -> std::string;

    /// Return the electric charge of the species.
    auto charge() const -> double;

    /// Return the molar mass of the species (in unit of kg/mol).
    auto molarMass() const -> double;

    /// Return the aggregate state of the species.
    auto aggregateState() const -> AggregateState;

    /// Return the elements that compose the species and their coefficients.
    auto elements() const -> const Elements&;

    /// Return the coefficient of an element in the species.
    auto elementCoefficient(const std::string& symbol) const -> double;

    /// Return the tags of the species (e.g., `organic`, `mineral`).
    auto tags() const -> const std::vector<std::string>&;

    /// Return the critical properties of the underlying substance of the species if available.
    auto criticalProps() const -> std::optional<SubstanceCriticalProps>;

    /// Return the attached data of the species whose type is known at runtime only.
    auto attachedData() const -> const std::any&;

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

// Custom specialization of std::hash for Reaktoro::Species
namespace std {

template<>
struct hash<Reaktoro::Species>
{
    std::size_t operator()(const Reaktoro::Species& s) const noexcept
    {
        return std::hash<std::string>{}(s.name());
    }
};

} // namespace std
