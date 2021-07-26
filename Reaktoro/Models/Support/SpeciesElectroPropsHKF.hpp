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
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

// Forward declarations
struct SpeciesElectroProps;
struct StandardThermoModelParamsHKF;
struct WaterThermoProps;

/// The *g* function state of in HKF model for computation of electrostatic properties of aqueous solutes.
struct gHKF
{
    /// The function *g* value.
    real g = {};

    /// The first-order partial derivative of function *g* with respect to temperature.
    real gT = {};

    /// The first-order partial derivative of function *g* with respect to pressure.
    real gP = {};

    /// The second-order partial derivative of function *g* with respect to temperature.
    real gTT = {};

    /// The second-order partial derivative of function *g* with respect to temperature and pressure.
    real gTP = {};

    /// The second-order partial derivative of function *g* with respect to pressure.
    real gPP = {};

    /// Compute the *g* function state at given temperature, pressure and water thermodynamic properties.
    static auto compute(real T, real P, const WaterThermoProps& wtp) -> gHKF;
};

/// Compute the electrostatic properties of an aqueous solute with given HKF *g* function state.
auto speciesElectroPropsHKF(const gHKF& gstate, const StandardThermoModelParamsHKF& params) -> SpeciesElectroProps;

} // namespace Reaktoro
