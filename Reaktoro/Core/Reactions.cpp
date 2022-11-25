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

#include "Reactions.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ReactionEquation.hpp>

namespace Reaktoro {

GeneralReaction::GeneralReaction()
{}

GeneralReaction::GeneralReaction(String const& equation)
{
    setName(equation);
    setEquation(equation);
}

auto GeneralReaction::setName(String const& name) -> GeneralReaction&
{
    reaction_name = name;
    return *this;
}

auto GeneralReaction::setEquation(String const& equation) -> GeneralReaction&
{
    reaction_equation = equation;
    return *this;
}

auto GeneralReaction::setRateModel(ReactionRateModel const& model) -> GeneralReaction&
{
    rate_model = model;
    rate_model_generator = {}; // clear any ReactionRateModelGenerator already set
    return *this;
}

auto GeneralReaction::setRateModel(ReactionRateModelGenerator const& model_generator) -> GeneralReaction&
{
    rate_model = {}; // clear any ReactionRateModel already set
    rate_model_generator = model_generator;
    return *this;
}

auto GeneralReaction::name() const -> String const&
{
    return reaction_name;
}

auto GeneralReaction::equation() const -> String const&
{
    return reaction_equation;
}

auto GeneralReaction::rateModel() const -> ReactionRateModel const&
{
    return rate_model;
}

auto GeneralReaction::rateModelGenerator() const -> ReactionRateModelGenerator  const&
{
    return rate_model_generator;
}

auto GeneralReaction::operator()(ReactionGeneratorArgs args) const -> Reaction
{
    // Ensure reaction name, equation, and rate model are given at the time of conversion
    errorif(reaction_name.empty(), "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction name. Use method GeneralReaction::setName to resolve this.");
    errorif(reaction_equation.empty(), "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction equation. Use method GeneralReaction::setEquation to resolve this.");
    errorif(!rate_model && !rate_model_generator, "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction rate model or a reaction rate model generator. Use method GeneralReaction::setRateModel to resolve this.");
    errorif(rate_model && rate_model_generator, "Converting a GeneralReaction object to a Reaction object requires either a non-empty reaction rate model or a reaction rate model generator. Do not use method GeneralReaction::setRateModel for both cases simultaneously.");

    // Construct the ReactionEquation object from string representation in `reaction_equation`
    ReactionEquation reaction_equation_obj(reaction_equation, args.species);

    // Collect the necessary data for reaction rate model generator.
    ReactionRateModelGeneratorArgs rargs{
        reaction_name,
        reaction_equation_obj,
        args.database,
        args.species,
        args.phases,
        args.surfaces,
    };

    // Resolve the reaction rate model of the reaction, either given or to be generated from a reaction rate model generator
    auto reaction_rate_model = rate_model ? rate_model : rate_model_generator(rargs);

    // Return the fully specified Reaction object
    return Reaction()
        .withName(reaction_name)
        .withEquation(reaction_equation_obj)
        .withRateModel(reaction_rate_model);
}

Reactions::Reactions()
{}

auto Reactions::convert(ReactionGeneratorArgs args) const -> Vec<Reaction>
{
    Vec<Reaction> reactions;

    for(auto const& fn : reaction_generators)
    {
        auto rxns = fn(args);
        reactions.insert(reactions.end(), rxns.begin(), rxns.end());
    }

    return reactions;
}

} // namespace Reaktoro
