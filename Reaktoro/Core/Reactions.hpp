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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionGenerator.hpp>

namespace Reaktoro {

/// The class used to represent a collection of reactions controlled kinetically.
/// @ingroup Core
class Reactions
{
public:
    /// Construct a Reactions object.
    Reactions();

    /// Construct a Reactions object with given generic reactions.
    /// @param reaction_generators The ReactionGenerator objects that will be converted into Reaction objects.
    template<typename... ReactionGenerators>
    Reactions(const ReactionGenerators&... reaction_generators)
    {
        static_assert(sizeof...(reaction_generators) > 0);
        addAux(reaction_generators...);
    }

    /// Add a reaction generator into the Reactions container.
    template<typename T>
    auto add(T const& item) -> void
    {
        static_assert(isBaseOf<ReactionGenerator, T> || isSame<T, Reaction>);

        if constexpr(isBaseOf<ReactionGenerator, T>)
        {
            rgeneratorfns.push_back([=](ChemicalSystem const& system) -> Vec<Reaction> {
                return item.convert(system); // item is a ReactionGenerator object
            });
        }
        else
        {
            rgeneratorfns.push_back([=](ChemicalSystem const& system) -> Vec<Reaction> {
                return { item }; // item is a Reaction object
            });
        }
    }

    /// Convert this Reactions object into a vector of Reaction objects.
    /// @param system The intermediate chemical system without attached reactions in which the reactions take place.
    auto convert(ChemicalSystem const& system) const -> Vec<Reaction>;

private:
    /// The ReactionGeneratorFn objects collected so far with each call to Reactions::add method.
    Vec<ReactionGeneratorFn> rgeneratorfns;

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
