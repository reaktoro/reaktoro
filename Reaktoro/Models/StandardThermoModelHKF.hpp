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

/// The parameters in the HKF model for calculating standard thermodynamic properties of aqueous solutes.
struct StandardThermoModelParamsHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in J/mol).
    Param Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in J/mol).
    Param Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in J/(mol·K)).
    Param Sr;

    /// The coefficient `a1` of the HKF equation of state of the aqueous solute (in J/(mol·Pa)).
    Param a1;

    /// The coefficient `a2` of the HKF equation of state of the aqueous solute (in J/mol).
    Param a2;

    /// The coefficient `a3` of the HKF equation of state of the aqueous solute (in (J·K)/(mol·Pa)).
    Param a3;

    /// The coefficient `a4` of the HKF equation of state of the aqueous solute (in (J·K)/mol).
    Param a4;

    /// The coefficient `c1` of the HKF equation of state of the aqueous solute (in J/(mol·K)).
    Param c1;

    /// The coefficient `c2` of the HKF equation of state of the aqueous solute (in (J·K)/mol).
    Param c2;

    /// The conventional Born coefficient of the aqueous solute at reference temperature 298.15 K and pressure 1 bar (in J/mol).
    Param wref;

    /// The electrical charge of the aqueous solute.
    real charge;

    /// The maximum temperature at which the HKF model can be applied for the substance (optional, in K).
    real Tmax;
};

/// Return a function that calculates thermodynamic properties of an aqueous solute using the HKF model.
auto StandardThermoModelHKF(const StandardThermoModelParamsHKF& params) -> StandardThermoModel;

} // namespace Reaktoro
