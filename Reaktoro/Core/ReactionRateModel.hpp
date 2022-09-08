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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Rate.hpp>
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

/// The type of functions that calculate reaction rates.
using ReactionRateModel = Model<Rate(ChemicalState const&)>;

/// The type of functions that construct a ReactionRateModel for a reaction.
/// @param reaction The reaction for which the reaction rate model is constructed.
/// @param phases The phases and their constituent species that form the chemical system.
using ReactionRateModelGenerator = Fn<ReactionRateModel(Reaction const& reaction, PhaseList const& phases)>;

} // namespace Reaktoro

//=========================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING CHEMICALPROPS
//=========================================================================

namespace Reaktoro {

template<typename T>
struct MemoizationTraits;

/// Specialize MemoizationTraits for ChemicalState.
template<>
struct MemoizationTraits<ChemicalState>
{
    using Type = ChemicalState;

    /// The type used instead to cache a ChemicalState object.
    using CacheType = Tuple<real, real, ArrayXr, ArrayXr>;

    static auto equal(const Tuple<real, real, ArrayXr, ArrayXr>& a, const ChemicalState& b)
    {
        auto const& [T, P, n, S] = a;
        return
            T == b.temperature() &&
            P == b.pressure() &&
            (n == b.speciesAmounts()).all() &&
            (S == b.surfaceAreas()).all();
    }

    static auto assign(Tuple<real, real, ArrayXr, ArrayXr>& a, const ChemicalState& b)
    {
        auto& [T, P, n, S] = a;
        T = b.temperature();
        P = b.pressure();
        n = b.speciesAmounts();
        S = b.surfaceAreas();
    }
};

} // namespace Reaktoro
