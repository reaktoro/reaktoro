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

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The parameters in the extended UNIQUAC model for calculating standard thermodynamic properties of substances.
struct StandardThermoModelParamsExtendedUNIQUAC
{
    /// The standard molar Gibbs energy of formation of the substance at reference temperature 298.15 K (in kJ/mol).
    real Gr;

    /// The standard molar enthalpy of formation of the substance at reference temperature 298.15 K (in kJ/mol).
    real Hr;

    /// The standard molar entropy of the substance at reference temperature 298.15 K (in J/(mol·K)).
    real Sr;

    /// The standard molar volume of the substance at reference temperature 298.15 K (in m³/mol).
    real Vr;

    /// The standard molar heat capacity (constant pressure) of the substance at reference temperature 298.15 K (in J/(mol·K)).
    real Cp;

    /// The parameter `a` in the standard heat capacity model (in J/(mol·K)).
    real a;

    /// The parameter `b` in the standard heat capacity model (in J/(mol·K²)).
    real b;

    /// The parameter `c` in the standard heat capacity model (in J/mol).
    real c;

    /// The parameter `α` in the standard heat capacity model (in 1/bar).
    real alpha;

    /// The parameter `β` in the standard heat capacity model (in 1/bar²).
    real beta;

    /// The parameter `Θ` in the standard heat capacity model (in K).
    real Theta = 200.0;
};

/// Return a function that calculates the thermodynamic properties of a substance using the extended UNIQUAC model.
auto StandardThermoModelExtendedUNIQUAC(StandardThermoModelParamsExtendedUNIQUAC const& params) -> StandardThermoModel;

} // namespace Reaktoro
