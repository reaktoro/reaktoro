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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionRateModel.hpp>

namespace Reaktoro {

class Database;
class PhaseList;
class SpeciesList;
class SurfaceList;

/// The data provided to a ReactionGenerator to construct a Reaction object.
/// @see ReactionGenerator, Reaction
/// @ingroup Core
struct ReactionGeneratorArgs
{
    /// The thermodynamic database used to construct the chemical system where the reaction belongs to.
    Database const& database;

    /// The species in the chemical system where the reaction belongs to.
    SpeciesList const& species;

    /// The phases in the chemical system where the reaction belongs to.
    PhaseList const& phases;

    /// The surfaces in the chemical system where the reaction belongs to.
    SurfaceList const& surfaces;
};

/// The function type for the generation of reactions with given species in the chemical system.
/// @param args The data provided to a ReactionGenerator to construct a Reaction object.
using ReactionGenerator = Fn<Vec<Reaction>(ReactionGeneratorArgs args)>;

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

    /// Set the reaction rate model generator of the reaction.
    /// Use this method to set a ReactionRateModelGenerator in case you need the
    /// ReactionRateModel of the reaction to be constructed later, when the
    /// chemical system is assembled.
    auto setRateModel(ReactionRateModelGenerator const& model_generator) -> GeneralReaction&;

    /// Set the reaction rate model generator of the reaction (equivalent to GeneralReaction::setRateModel).
    auto set(ReactionRateModelGenerator const& model_generator) -> GeneralReaction&;

    /// Return the name of the reaction.
    auto name() const -> String const&;

    /// Return the reaction equation of the reaction.
    auto equation() const -> String const&;

    /// Return the reaction rate model of the reaction.
    auto rateModel() const -> ReactionRateModel const&;

    /// Return the reaction rate model generator of the reaction.
    auto rateModelGenerator() const -> ReactionRateModelGenerator const&;

    /// Convert this GeneralReaction object into a Reaction object.
    auto operator()(ReactionGeneratorArgs args) const -> Reaction;

private:
    /// The name of the reaction.
    String reaction_name;

    /// The reaction equation of the reaction as a formatted string.
    String reaction_equation;

    /// The rate model of the reaction.
    ReactionRateModel rate_model;

    /// The rate model generator of the reaction.
    ReactionRateModelGenerator rate_model_generator;
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
            reaction_generators.push_back([=](ReactionGeneratorArgs args) -> Vec<Reaction> { return { item }; }); // item is a Reaction object
        }
        else if constexpr(isConvertible<T, GeneralReaction>)
        {
            reaction_generators.push_back([=](ReactionGeneratorArgs args) -> Vec<Reaction> { return { item(args) }; }); // item is a GeneralReaction object; use operator()(ReactionGeneratorArgs) to convert to Reaction
        }
        else
        {
            reaction_generators.push_back(item); // item is already a ReactionGenerator
            errorif(!reaction_generators.back(), "Expecting an object in Reactions::add call that produces an initialized ReactionGenerator function object.");
        }
    }

    /// Convert this Reactions object into a vector of Reaction objects.
    /// @param args The data provided to each ReactionGenerator to construct the Reaction objects.
    auto convert(ReactionGeneratorArgs args) const -> Vec<Reaction>;

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
