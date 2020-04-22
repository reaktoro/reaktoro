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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// The type used to store HKF parameters based on SUPCRT databases for the aqueous solvent (water).
struct SupcrtParamsAqueousSolventHKF
{
};

/// The type used to store HKF parameters based on SUPCRT databases for an aqueous solute.
struct SupcrtParamsAqueousSoluteHKF
{
    /// The name of the aqueous solute.
    String name;

    /// The electrical charge of the aqueous solute
    real charge = {};

    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    real Gf = {};

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    real Hf = {};

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    real Sr = {};

    /// The coefficient a1 of the HKF equation of state of the aqueous solute (in units of cal/(mol*bar))
    real a1 = {};

    /// The coefficient a2 of the HKF equation of state of the aqueous solute (in units of cal/mol)
    real a2 = {};

    /// The coefficient a3 of the HKF equation of state of the aqueous solute (in units of (cal*K)/(mol*bar))
    real a3 = {};

    /// The coefficient a4 of the HKF equation of state of the aqueous solute (in units of (cal*K)/mol)
    real a4 = {};

    /// The coefficient c1 of the HKF equation of state of the aqueous solute (in units of cal/(mol*K))
    real c1 = {};

    /// The coefficient c2 of the HKF equation of state of the aqueous solute (in units of (cal*K)/mol)
    real c2 = {};

    /// The conventional Born coefficient of the aqueous solute at reference temperature 298.15 K and pressure 1 bar (in units of cal/mol)
    real wref = {};
};

/// The type used to store Maier-Kelly parameters based on SUPCRT databases for gaseous/liquid species.
struct SupcrtParamsMaierKelly
{
    /// The name of the species.
    String name;

    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    real Gf = {};

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    real Hf = {};

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    real Sr = {};

    /// The coefficient a of the Maier-Kelly equation of state of the gaseous/liquid species (in units of cal/(mol*K))
    real a = {};

    /// The coefficient b of the Maier-Kelly equation of state of the gaseous/liquid species (in units of cal/(mol*K^2))
    real b = {};

    /// The coefficient c of the Maier-Kelly equation of state of the gaseous/liquid species (in units of (cal*K)/mol)
    real c = {};

    /// The maximum temperature at which the Maier-Kelly equation of state can be applied for the gaseous species (in units of K)
    real Tmax = {};
};

/// The type used to store Maier-Kelly-HKF parameters based on SUPCRT databases for mineral species.
struct SupcrtParamsMaierKellyHKF
{
    /// The name of the species.
    String name;

    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    real Gf = {};

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    real Hf = {};

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    real Sr = {};

    /// The standard molal volume of the mineral species at reference temperature and pressure (in units of cm^3/mol)
    real Vr = {};

    /// The number of phase transitions that the mineral may undergo
    Index nptrans = 0;

    /// The coefficients ai of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of cal/(mol*K))
    Vec<real> a;

    /// The coefficients bi of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of cal/(mol*K^2))
    Vec<real> b;

    /// The coefficients ci of the Maier-Kelly-HKF equation of state of the mineral species for each phase region (in units of (cal*K)/mol)
    Vec<real> c;

    /// The temperatures at which the mineral experiences phase transition along the line of reference pressure (in units of K)
    Vec<real> Ttr;

    /// The change in the standard molal enthalpy of each mineral phase transition (in units of cal/mol)
    Vec<real> Htr;

    /// The change in the standard molal volume of each mineral phase transition (in units of cm^3/mol)
    Vec<real> Vtr;

    /// The Clapeyron slote at each mineral phase transition (in units of bar/K)
    Vec<real> dPdTtr;

    /// The maximum temperature at which the Maier-Kelly-HKF equation of state can be applied for the mineral (in units of K)
    real Tmax = {};
};

} // namespace Reaktoro
