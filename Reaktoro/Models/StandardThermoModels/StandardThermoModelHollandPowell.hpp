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

/// The parameters in the Holland-Powell model for calculating standard thermodynamic properties of fluid and mineral species.
struct StandardThermoModelParamsHollandPowell
{
    /// The apparent standard molal Gibbs free energy of formation of the substance from its elements (in J/mol).
    real Gf;

    /// The apparent standard molal enthalpy of formation of the substance from its elements (in J/mol).
    real Hf;

    /// The standard molal entropy of the substance at reference temperature and pressure (in J/(mol·K)).
    real Sr;

    /// The standard molal volume of the substance at reference temperature and pressure (in m³/mol).
    real Vr;

    /// The coefficient `a` of the Holland and Powell (2011) thermodynamic model (in J/(mol·K)).
    real a;

    /// The coefficient `b` of the Holland and Powell (2011) thermodynamic model (in J/(mol·K²)).
    real b;

    /// The coefficient `c` of the Holland and Powell (2011) thermodynamic model (in (J·K)/mol).
    real c;

    /// The coefficient `d` of the Holland and Powell (2011) thermodynamic model (in J/(mol·K½)).
    real d;

    /// The coefficient `α0` of the Holland and Powell (2011) thermodynamic model (in 1/K).
    real alpha0;

    /// The coefficient `κ0` of the Holland and Powell (2011) thermodynamic model (in Pa).
    real kappa0;

    /// The coefficient `κ0'` of the Holland and Powell (2011) thermodynamic model (dimensionless).
    real kappa0p;

    /// The coefficient `κ0''` of the Holland and Powell (2011) thermodynamic model (in 1/Pa).
    real kappa0pp;

    /// The number of atoms in the chemical formula of the mineral.
    real numatoms;

    /// The maximum temperature at which the Holland and Powell (2011) model can be applied for the substance (optional, in K).
    real Tmax;

    // TODO: The following data, although provided in SUPCRTBL for some
    // minerals, is from HP98, and not HP11, and thus not actually used. Decide
    // to remove this.

    // /// The critical temperature of the mineral at 1 bar (in K).
    // real Tcr;

    // /// The entropy of disordering of the mineral at the critical temperature above (in J/(mol·K)).
    // real Smax;

    // /// The volume of disordering of the mineral at the critical temperature above (in m³/mol).
    // real Vmax;
};

/// Return a function that calculates thermodynamic properties of a fluid or mineral species using the Holland-Powell model.
auto StandardThermoModelHollandPowell(const StandardThermoModelParamsHollandPowell& params) -> StandardThermoModel;

} // namespace Reaktoro
