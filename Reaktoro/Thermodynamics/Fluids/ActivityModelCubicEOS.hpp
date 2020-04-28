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

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Thermodynamics/Fluids/CubicEOS.hpp>

namespace Reaktoro {

// Forward declarations
class GeneralMixture;

/// The options for the activity model based on cubic equation of state.
struct ActivityModelOptionsCubicEOS
{
    /// The fluid type for which the equation of state should be confifured.
    CubicEOSFluidType fluidtype = CubicEOSFluidType::Vapor;

    /// The cubic equation of state model to be used.
    CubicEOSModel model = CubicEOSModel::PengRobinson;

    /// The function that calculates interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
    CubicEOSInteractionParamsFn interaction_params_fn;

    /// The method to identify whether liquid or vapor phases is stable.
    PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;
};

/// Return an activity model for a fluid, liquid or gaseous, based on a cubic equation of state.
auto activityModelCubicEOS(const GeneralMixture& mixture, ActivityModelOptionsCubicEOS options) -> ActivityPropsFn;

} // namespace Reaktoro
