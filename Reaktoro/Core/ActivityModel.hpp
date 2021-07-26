// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The arguments in an function for calculation of activity properties of a phase.
/// @see ActivityModel, ActivityProps
struct ActivityArgs
{
    /// The temperature for the calculation (in K).
    const real& T;

    /// The pressure for the calculation (in Pa).
    const real& P;

    /// The mole fractions of the species in the phase.
    ArrayXrConstRef x;
};

// Declare this so that Model undestand ActivityPropsRef as reference type for ActivityProps instead of ActivityProps&.
REAKTORO_DEFINE_REFERENCE_TYPE_OF(ActivityProps, ActivityPropsRef);

/// The function type for the calculation of activity and excess thermodynamic properties of a phase.
using ActivityModel = Model<ActivityProps(ActivityArgs)>;

/// The type for functions that construct an ActivityModel for a phase.
/// @param species The species in the phase.
using ActivityModelGenerator = Fn<ActivityModel(const SpeciesList& species)>;

/// Return an activity model resulting from chaining other activity models.
auto chain(const Vec<ActivityModelGenerator>& models) -> ActivityModelGenerator;

/// Return an activity model resulting from chaining other activity models.
auto chain(const ActivityModelGenerator& model) -> ActivityModelGenerator;

/// Return an activity model resulting from chaining other activity models.
template<typename... Models>
auto chain(const ActivityModelGenerator& model, const Models&... models) -> ActivityModelGenerator
{
    Vec<ActivityModelGenerator> vec = {model, models...};
    return chain(vec);
}

} // namespace Reaktoro
