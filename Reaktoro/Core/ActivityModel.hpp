// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ActivityProps.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The base class for all other activity model classes.
class ActivityModel
{
public:
    /// Build the function for activity and thermodynamic excesss property calculations of a phase.
    /// @param species The species in the phase.
    virtual auto build(const SpeciesList& species) const -> ActivityPropsFn = 0;

    /// Convert the derived activity model object into an activity property function.
    auto operator()(const SpeciesList& species) -> ActivityPropsFn { return build(species); }
};

/// The function type for the creation of activity model function of a phase.
using ActivityModelFn = Fn<ActivityPropsFn(const SpeciesList&)>;

} // namespace Reaktoro
