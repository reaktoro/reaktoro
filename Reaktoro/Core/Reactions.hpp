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
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

/// The function type for the generation of reactions with given phases and their species in the chemical system.
using ReactionGenerator = Fn<Vec<Reaction>(PhaseList const&)>;

/// The class used to represent a collection of reactions controlled kinetically.
/// @ingroup Core
class Reactions
{
public:
    /// Construct a Reactions object.
    Reactions();

    /// Construct a Reactions object with given reaction generators.
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
        static_assert(isConvertible<T, ReactionGenerator> || isConvertible<T, Reaction>);

        if constexpr(isConvertible<T, ReactionGenerator>)
        {
            rgenerators.push_back([=](PhaseList const& phases) -> Vec<Reaction> {
                return item(phases); // item is a ReactionGenerator object
            });
        }
        else
        {
            rgenerators.push_back([=](PhaseList const& phases) -> Vec<Reaction> {
                return { Reaction(item) }; // item is convertible to Reaction object
            });
        }
    }

    /// Convert this Reactions object into a vector of Reaction objects.
    /// @param phases The phases and their species composing the chemical system where the reactions will take place.
    auto convert(PhaseList const& phases) const -> Vec<Reaction>;

private:
    /// The ReactionGenerator objects collected so far with each call to Reactions::add method.
    Vec<ReactionGenerator> rgenerators;

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
