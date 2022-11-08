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
    return *this;
}

auto GeneralReaction::set(ReactionRateModel const& model) -> GeneralReaction&
{
    return setRateModel(model);
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

auto GeneralReaction::operator()(SpeciesList const& species) const -> Reaction
{
    errorif(reaction_name.empty(), "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction name. Use method GeneralReaction::setName to resolve this.");
    errorif(reaction_equation.empty(), "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction equation. Use method GeneralReaction::setEquation to resolve this.");
    errorif(!rate_model.initialized(), "Converting a GeneralReaction object to a Reaction object requires a non-empty reaction rate model. Use method GeneralReaction::setRateModel to resolve this.");

    return Reaction()
        .withName(reaction_name)
        .withEquation(ReactionEquation(reaction_equation, species))
        .withRateModel(rate_model);
}

Reactions::Reactions()
{}

auto Reactions::convert(SpeciesList const& species) const -> Vec<Reaction>
{
    Vec<Reaction> reactions;

    for(auto const& fn : reaction_generators)
    {
        auto rxns = fn(species);
        reactions.insert(reactions.end(), rxns.begin(), rxns.end());
    }

    return reactions;
}

} // namespace Reaktoro
