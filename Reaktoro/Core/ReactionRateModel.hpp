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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/Rate.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

/// The type of functions that calculate reaction rates.
using ReactionRateModel = Model<Rate(ChemicalProps const&)>;

/// The type of functions that construct a ReactionRateModel for a reaction.
/// @param reaction The reaction for which the reaction rate model is constructed.
/// @param phases The phases and their constituent species that form the chemical system.
using ReactionRateModelGenerator = Fn<ReactionRateModel(Reaction const& reaction, Phases const& phases)>;

} // namespace Reaktoro

//=========================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING CHEMICALPROPS
//=========================================================================

namespace Reaktoro {

template<typename T>
struct MemoizationTraits;

/// Specialize MemoizationTraits for ChemicalProps.
template<>
struct MemoizationTraits<ChemicalProps>
{
    using Type = ChemicalProps;

    /// The type used instead to cache a ChemicalProps object.
    using CacheType = Tuple<real, real, ArrayXr>;

    static auto equal(const Tuple<real, real, ArrayXr>& a, const ChemicalProps& b)
    {
        auto const& [T, P, n] = a;
        return T == b.temperature() && P == b.pressure() && (n == b.speciesAmounts()).all();
    }

    static auto assign(Tuple<real, real, ArrayXr>& a, const ChemicalProps& b)
    {
        auto& [T, P, n] = a;
        T = b.temperature();
        P = b.pressure();
        n = b.speciesAmounts();
    }
};

} // namespace Reaktoro
