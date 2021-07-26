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

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Thermodynamics/Fluids/CubicEOS.hpp>

namespace Reaktoro {

/// The parameters for the activity model based on cubic equations of state.
struct ActivityModelCubicEOSParams
{
    /// The function that calculates interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
    CubicEOSInteractionParamsFn interaction_params_fn = nullptr;

    /// The method to identify whether liquid or vapor phases is stable.
    PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;
};

/// Return the activity model for fluid phases based on the Van der Waals cubic equation of state.
auto ActivityModelVanDerWaals(ActivityModelCubicEOSParams params = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Redlich-Kwong cubic equation of state.
auto ActivityModelRedlichKwong(ActivityModelCubicEOSParams params = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Soave-Redlich-Kwong cubic equation of state.
auto ActivityModelSoaveRedlichKwong(ActivityModelCubicEOSParams params = {}) -> ActivityModelGenerator;

/// Return the activity model for fluid phases based on the Peng-Robinson cubic equation of state.
auto ActivityModelPengRobinson(ActivityModelCubicEOSParams params = {}) -> ActivityModelGenerator;

} // namespace Reaktoro
