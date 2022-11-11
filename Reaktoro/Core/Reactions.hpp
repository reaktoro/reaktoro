// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The function type for the generation of reactions with given species in the chemical system.
/// @param species The species composing the chemical system where the reactions will take place.
using ReactionGenerator = Fn<Vec<Reaction>(SpeciesList const& species)>;

/// Used to define a general reaction.
/// @ingroup Core
class GeneralReaction
{
public:
    /// Construct a default GeneralReaction object.
    GeneralReaction();

    /// Construct a GeneralReaction object with given reaction equation as a formatted string.
    /// @copydetails GeneralReaction::setEquation
    explicit GeneralReaction(String const& equation);

    /// Set the unique name of the reaction.
    auto setName(String const& name) -> GeneralReaction&;

    /// Set the equation of the reaction as a formatted string.
    /// Below are examples of formatted strings representing reaction equations:
    /// ~~~
    /// GeneralReaction("CO2(g) + H2O = H+ + HCO3-");
    /// GeneralReaction("Calcite + H+ = Ca++ + HCO3-");
    /// GeneralReaction("Dolomite + 2*H+ = Ca++ + Mg++ + 2*HCO3-");
    /// GeneralReaction("Quartz");
    /// ~~~
    /// Note that unity stoichiometric coefficients can be ommited from the equation. The operator
    /// `*` must be used when this is not the case. Also note that equations with a single reactant
    /// name is possible, as in the Quartz reaction above. This exists to permit reaction rates to
    /// be defined for a single species instead of a normal stoichiometrically balanced reaction.
    auto setEquation(String const& equation) -> GeneralReaction&;

    /// Set the reaction rate model of the reaction.
    auto setRateModel(ReactionRateModel const& model) -> GeneralReaction&;

    /// Set the reaction rate model of the reaction (equivalent to GeneralReaction::setRateModel).
    auto set(ReactionRateModel const& model) -> GeneralReaction&;

    /// Return the name of the reaction.
    auto name() const -> String const&;

    /// Return the reaction equation of the reaction.
    auto equation() const -> String const&;

    /// Return the reaction rate model of the reaction.
    auto rateModel() const -> ReactionRateModel const&;

    /// Convert this GeneralReaction object into a Reaction object.
    auto operator()(SpeciesList const& species) const -> Reaction;

private:
    /// The name of the reaction.
    String reaction_name;

    /// The reaction equation of the reaction as a formatted string.
    String reaction_equation;

    /// The rate model of the reaction.
    ReactionRateModel rate_model;
};

/// Used to represent a collection of reactions controlled kinetically.
/// @ingroup Core
class Reactions
{
public:
    /// Construct a Reactions object.
    Reactions();

    /// Construct a Reactions object with given Reaction, GeneralReaction, or ReactionGenerator objects.
    /// @param reactions The objects of type Reaction, GeneralReaction, or ReactionGenerator.
    template<typename... ReactionConvertible>
    explicit Reactions(ReactionConvertible const&... reactions)
    {
        static_assert(sizeof...(reactions) > 0);
        addAux(reactions...);
    }

    /// Add a reaction generator into the Reactions container.
    template<typename T>
    auto add(T const& item) -> void
    {
        static_assert(
            isConvertible<T, Reaction> ||
            isConvertible<T, GeneralReaction> ||
            isConvertible<T, ReactionGenerator>);

        if constexpr(isConvertible<T, Reaction>)
        {
            reaction_generators.push_back([=](SpeciesList const& species) -> Vec<Reaction> { return { item }; }); // item is a Reaction object
        }
        else if constexpr(isConvertible<T, GeneralReaction>)
        {
            reaction_generators.push_back([=](SpeciesList const& species) -> Vec<Reaction> { return { item(species) }; }); // item is a GeneralReaction object; use operator()(SpeciesList) to convert to Reaction
        }
        else
        {
            reaction_generators.push_back(item); // item is already a ReactionGenerator
        }
    }

    /// Convert this Reactions object into a vector of Reaction objects.
    /// @param species The species composing the chemical system where the reactions will take place.
    auto convert(SpeciesList const& species) const -> Vec<Reaction>;

private:
    /// The ReactionGenerator objects collected so far with each call to Reactions::add method.
    Vec<ReactionGenerator> reaction_generators;

    /// Add one or more ReactionGenerator or Reaction objects into the Reactions container.
    template<typename Arg, typename... Args>
    auto addAux(const Arg& arg, const Args&... args) -> void
    {
        add(arg);
        if constexpr (sizeof...(Args) > 0)
            addAux(args...);
    }
};

} // namespace Reaktoro
