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
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

/// The parameters in the Maier-Kelley-HKF model for calculating standard thermodynamic properties of mineral species.
struct StandardThermoModelParamsMineralHKF
{
    /// The apparent standard Gibbs free energy of formation of the mineral substance from its elements (in J/mol).
    Param Gf;

    /// The apparent standard enthalpy of formation of the mineral substance from its elements (in J/mol).
    Param Hf;

    /// The standard entropy of the mineral substance at reference temperature and pressure (in J/(mol·K)).
    Param Sr;

    /// The standard volume of the mineral substance at reference temperature and pressure (in m³/mol).
    Param Vr;

    /// The number of phase transitions considered in the mineral.
    Index ntr;

    /// The coefficients `a(i)` of the Maier-Kelley-HKF mineral thermodynamic model for each phase region (in J/(mol·K)).
    Vec<Param> a;

    /// The coefficients `b(i)` of the Maier-Kelley-HKF mineral thermodynamic model for each phase region (in J/(mol·K²)).
    Vec<Param> b;

    /// The coefficients `c(i)` of the Maier-Kelley-HKF mineral thermodynamic model for each phase region (in (J·K)/mol).
    Vec<Param> c;

    /// The temperatures at which the mineral experiences phase transition along the line of reference pressure (in K).
    Vec<Param> Ttr;

    /// The change in the standard enthalpy of each mineral phase transition (in J/mol).
    Vec<Param> Htr;

    /// The change in the standard volume of each mineral phase transition (in m³/mol).
    Vec<Param> Vtr;

    /// The Clapeyron slote at each mineral phase transition (in Pa/K).
    Vec<Param> dPdTtr;

    /// The maximum temperature at which the Maier-Kelley-HKF model can be applied for the substance (optional, in K).
    real Tmax;
};

/// Return a function that calculates thermodynamic properties of a species using the Maier-Kelley-HKF model.
auto StandardThermoModelMineralHKF(const StandardThermoModelParamsMineralHKF& params) -> StandardThermoModel;

} // namespace Reaktoro
