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

// C++ includes
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ActivityProps.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The arguments in an function for calculation of activity properties of a phase.
/// @see ActivityModel, ActivityProps
struct ActivityModelArgs
{
    /// The temperature for the calculation (in K).
    real const& T;

    /// The pressure for the calculation (in Pa).
    real const& P;

    /// The mole fractions of the species in the phase.
    ArrayXrConstRef x;
};

// Declare this so that Model understands ActivityPropsRef as reference type for ActivityProps instead of ActivityProps&.
REAKTORO_DEFINE_REFERENCE_TYPE_OF(ActivityProps, ActivityPropsRef);

/// The function type for the calculation of activity and corrective thermodynamic properties of a phase.
using ActivityModel = Model<ActivityProps(ActivityModelArgs)>;

/// The type for functions that construct an ActivityModel for a phase.
/// @param species The species in the phase.
using ActivityModelGenerator = Fn<ActivityModel(SpeciesList const& species)>;

/// Return an activity model resulting from chaining other activity models.
auto chain(const Vec<ActivityModelGenerator>& models) -> ActivityModelGenerator;

/// Return an activity model resulting from chaining other activity models.
auto chain(ActivityModelGenerator const& model) -> ActivityModelGenerator;

/// Return an activity model resulting from chaining other activity models.
template<typename... Models>
auto chain(ActivityModelGenerator const& model, Models const&... models) -> ActivityModelGenerator
{
    Vec<ActivityModelGenerator> vec = {model, models...};
    return chain(vec);
}

} // namespace Reaktoro

//=========================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING ACTIVITYMODELARGS
//=========================================================================

namespace Reaktoro {

template<typename T>
struct MemoizationTraits;

/// Specialize MemoizationTraits for ActivityModelArgs.
template<>
struct MemoizationTraits<ActivityModelArgs>
{
    using Type = ActivityModelArgs;

    /// The type used instead to cache an ActivityModelArgs object.
    using CacheType = Tuple<real, real, ArrayXr>;

    static auto equal(Tuple<real, real, ArrayXr> const& a, ActivityModelArgs const& b)
    {
        auto const& [T, P, x] = a;
        return T == b.T && P == b.P && (x == b.x).all();
    }

    static auto assign(Tuple<real, real, ArrayXr>& a, ActivityModelArgs const& b)
    {
        auto& [T, P, x] = a;
        T = b.T;
        P = b.P;
        x = b.x;
    }
};

} // namespace Reaktoro
