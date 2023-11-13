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

/// The parameters in the extended UNIQUAC model for calculating standard thermodynamic properties of substances.
struct StandardThermoModelParamsExtendedUNIQUAC
{
    /// The standard molar Gibbs energy of formation of the substance at reference temperature 298.15 K (in kJ/mol).
    Param Gr;

    /// The standard molar enthalpy of formation of the substance at reference temperature 298.15 K (in kJ/mol).
    Param Hr;

    /// The standard molar entropy of the substance at reference temperature 298.15 K (in J/(mol·K)).
    Param Sr;

    /// The standard molar volume of the substance at reference temperature 298.15 K (in m³/mol).
    Param Vr;

    /// The standard molar heat capacity (constant pressure) of the substance at reference temperature 298.15 K (in J/(mol·K)).
    Param Cp;

    /// The parameter `a` in the standard heat capacity model (in J/(mol·K)).
    Param a;

    /// The parameter `b` in the standard heat capacity model (in J/(mol·K²)).
    Param b;

    /// The parameter `c` in the standard heat capacity model (in J/mol).
    Param c;

    /// The parameter `α` in the standard heat capacity model (in 1/bar).
    Param alpha;

    /// The parameter `β` in the standard heat capacity model (in 1/bar²).
    Param beta;

    /// The parameter `Θ` in the standard heat capacity model (in K).
    Param Theta = 200.0;
};

/// Return a function that calculates the thermodynamic properties of a substance using the extended UNIQUAC model.
auto StandardThermoModelExtendedUNIQUAC(StandardThermoModelParamsExtendedUNIQUAC const& params) -> StandardThermoModel;

} // namespace Reaktoro
