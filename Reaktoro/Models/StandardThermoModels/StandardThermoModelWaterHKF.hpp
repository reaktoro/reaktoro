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
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The parameters in the HKF model for calculating standard thermodynamic properties of the aqueous solvent (water).
/// The default values below were taken from Helgeson and Kirkham (1974), page 1098.
struct StandardThermoModelParamsWaterHKF
{
    /// The temperature of liquid water at the triple point (in K).
    real Ttr = 273.16;

    /// The molar entropy of liquid water at the triple point (in J/(mol·K)).
    real Str = 63.312288;

    /// The molar Gibbs energy of liquid water at the triple point (in J/mol).
    real Gtr = -235517.36;

    /// The molar enthalpy of liquid water at the triple point (in J/mol).
    real Htr = -287721.128;
};

/// Return a function that calculates thermodynamic properties of the aqueous solvent (water) using the HKF model.
auto StandardThermoModelWaterHKF(const StandardThermoModelParamsWaterHKF& params) -> StandardThermoModel;

} // namespace Reaktoro
