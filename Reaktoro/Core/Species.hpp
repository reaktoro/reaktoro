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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ElementalComposition.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

/// A type used to represent a chemical species and its attributes.
class Species
{
public:
    /// Construct a default Species object.
    Species();

    /// Construct a Species object with given chemical formula (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    explicit Species(String formula);

    /// Return a deep copy of this Species object.
    auto clone() const -> Species;

    /// Return a duplicate of this Species object with new name that uniquely identifies this species.
    auto withName(String name) const -> Species;

    /// Return a duplicate of this Species object with new formula (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    auto withFormula(String formula) const -> Species;

    /// Return a duplicate of this Species object with new case insensitive substance name (e.g. `WATER`, `CARBON-MONOXIDE`).
    auto withSubstance(String substance) const -> Species;

    /// Return a duplicate of this Species object with new elemental composition.
    auto withElements(ElementalComposition elements) const -> Species;

    /// Return a duplicate of this Species object with new electric charge attribute.
    auto withCharge(double charge) const -> Species;

    /// Return a duplicate of this Species object with new aggregate state.
    auto withAggregateState(AggregateState option) const -> Species;

    /// Return a duplicate of this Species object with new formation reaction and new standard thermodynamic model.
    /// This method will also set the standard thermodynamic model of the species
    /// using the standard thermodynamic model assigned to the formation reaction.
    /// Use this method thus to assign a standard thermodynamic model to the
    /// Species object instead of using method @ref withStandardThermoPropsFn.
    auto withFormationReaction(FormationReaction reaction) const -> Species;

    /// Return a duplicate of this Species object with new standard thermodynamic model.
    /// This method exists for convenience only. Its use results in a standard
    /// thermodynamic model for this species in which its standard Gibbs energy
    /// is constant. All other standard thermodynamic properties are set to
    /// zero. For a more complete standard thermodynamic model, use method @ref
    /// withStandardThermoModel or @ref withFormationReaction, in case the
    /// thermodynamic model is based on reaction properties.
    /// @param G0 The constant standard Gibbs energy of the species (in J/mol).
    auto withStandardGibbsEnergy(real G0) const -> Species;

    /// Return a duplicate of this Species object with new standard thermodynamic model.
    /// This method assigns a standard thermodynamic model for the computation
    /// of standard thermodynamic properties of the species at given
    /// temperature and pressure. Alternatively, methods @ref
    /// withStandardGibbsEnergy and @ref withFormationReaction can be used to
    /// indirectly assign a standard thermodynamic model to this species.
    auto withStandardThermoPropsFn(const StandardThermoPropsFn& fn) const -> Species;

    /// Return a duplicate of this Species object with new tags attribute.
    auto withTags(const Strings& tags) const -> Species;

    /// Return a duplicate of this Species object with new attached data whose type is known at runtime only.
    auto withAttachedData(Any data) const -> Species;

    /// Return the name that uniquely identifies this species if provided, otherwise, its formula.
    auto name() const -> String;

    /// Return the chemical formula of the species.
    auto formula() const -> ChemicalFormula;

    /// Return the name of the underlying substance of the species if provided, otherwise, its formula without any suffix.
    auto substance() const -> String;

    /// Return the elements that compose the species and their coefficients.
    auto elements() const -> const ElementalComposition&;

    /// Return the electric charge of the species.
    auto charge() const -> double;

    /// Return the aggregate state of the species.
    auto aggregateState() const -> AggregateState;

    /// Return the formation reaction of the species.
    auto reaction() const -> const FormationReaction&;

    /// Return the function that computes the standard thermodynamic properties of the species.
    auto standardThermoPropsFn() const -> const StandardThermoPropsFn&;

    /// Return the tags of the species (e.g., `organic`, `mineral`).
    auto tags() const -> const Strings&;

    /// Return the attached data of the species whose type is known at runtime only.
    auto attachedData() const -> const Any&;

    /// Return the molar mass of the species (in unit of kg/mol).
    auto molarMass() const -> double;

    /// Return the standard thermodynamic properties of the species at given temperature (in K) and pressure (in Pa).
    auto props(real T, real P) const -> StandardThermoProps;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Return true if a Species object is less than another for sorting reasons.
auto operator<(const Species& lhs, const Species& rhs) -> bool;

/// Return true if two Species objects have the same name.
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
