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

/// The parameters in the Maier-Kelley model for calculating standard thermodynamic properties of fluid and solid species.
struct StandardThermoModelParamsMaierKelley
{
    /// The apparent standard molar Gibbs free energy of formation of the substance from its elements (in J/mol).
    real Gf;

    /// The apparent standard molar enthalpy of formation of the substance from its elements (in J/mol).
    real Hf;

    /// The standard molar entropy of the substance at reference temperature and pressure (in J/(mol·K)).
    real Sr;

    /// The standard molar volume of the mineral substance at reference temperature and pressure (in unit of m³/mol)
    real Vr;

    /// The coefficient `a` of the Maier-Kelley model (in J/(mol·K)).
    real a;

    /// The coefficient `b` of the Maier-Kelley model (in J/(mol·K²)).
    real b;

    /// The coefficient `c` of the Maier-Kelley model (in (J·K)/mol).
    real c;

    /// The maximum temperature at which the Maier-Kelley model can be applied for the substance (optional, in K).
    real Tmax;
};

/// Return a function that calculates thermodynamic properties of a species using the Maier-Kelley model.
auto StandardThermoModelMaierKelley(const StandardThermoModelParamsMaierKelley& params) -> StandardThermoModel;

} // namespace Reaktoro
