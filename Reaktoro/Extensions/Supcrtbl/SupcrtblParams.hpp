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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// The SUPCRTBL parameters for calculating standard thermodynamic properties of aqueous solvent water using HKF model.
/// Default values taken from Helgeson and Kirkham (1974), page 1098.
struct SupcrtblParamsAqueousSolventHKF
{
    /// The temperature of liquid water at the triple point (in K).
    real Ttr = 273.16;

    /// The molar entropy of liquid water at the triple point (in J/(mol·K)).
    real Str = 63.312288;

    /// The molar Gibbs energy of liquid water at the triple point (in kJ/mol).
    real Gtr = -235.51736;

    /// The molar enthalpy of liquid water at the triple point (in kJ/mol).
    real Htr = -287.721128;

    /// The molar internal energy of liquid water at the triple point (in kJ/mol).
    real Utr = -284.039208;

    /// The molar Helmholtz energy of liquid water at the triple point (in kJ/mol).
    real Atr = -231.85636;
};

/// The SUPCRTBL parameters for calculating standard thermodynamic properties of aqueous solutes using HKF model.
struct SupcrtblParamsAqueousSoluteHKF
{
    /// The name of the aqueous solute.
    String name;

    /// The electrical charge of the aqueous solute
    real charge = {};

    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in kJ/mol)
    real Gf = {};

    /// The apparent standard molal enthalpy of formation of the species from its elements (in kJ/mol)
    real Hf = {};

    /// The standard molal entropy of the species at reference temperature and pressure (in J/(mol·K))
    real Sr = {};

    /// The coefficient `a1` of the HKF equation of state of the aqueous solute (in kJ/(mol·bar)) (multiplied by 1e+1)
    real a1 = {};

    /// The coefficient `a2` of the HKF equation of state of the aqueous solute (in kJ/mol) (multiplied by 1e-2)
    real a2 = {};

    /// The coefficient `a3` of the HKF equation of state of the aqueous solute (in (kJ·K)/(mol·bar))
    real a3 = {};

    /// The coefficient `a4` of the HKF equation of state of the aqueous solute (in (kJ·K)/mol) (multiplied by 1e-4)
    real a4 = {};

    /// The coefficient `c1` of the HKF equation of state of the aqueous solute (in kJ/(mol·K))
    real c1 = {};

    /// The coefficient `c2` of the HKF equation of state of the aqueous solute (in (kJ·K)/mol) (multiplied by 1e-4)
    real c2 = {};

    /// The conventional Born coefficient of the aqueous solute at reference temperature 298.15 K and pressure 1 bar (in kJ/mol) (multiplied by 1e-5)
    real wref = {};
};

/// The SUPCRTBL parameters for calculating standard thermodynamic properties of gases/liquids using Holland and Powell (2011) model.
struct SupcrtblParamsFluidHollandPowell
{
    /// The apparent standard molal Gibbs free energy of formation of the substance from its elements (in kJ/mol)
    real Gf;

    /// The apparent standard molal enthalpy of formation of the substance from its elements (in kJ/mol)
    real Hf;

    /// The standard molal entropy of the substance at reference temperature and pressure (in J/(mol·K))
    real Sr;

    /// The standard molal volume of the substance at reference temperature and pressure (in J/(bar·mol))
    real Vr;

    /// The coefficient `a` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K))
    real a;

    /// The coefficient `b` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K²)) (multiplied by 1e+5)
    real b;

    /// The coefficient `c` of the Holland and Powell (2011) thermodynamic model (in (kJ·K)/mol)
    real c;

    /// The coefficient `d` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K½))
    real d;

    /// The maximum temperature at which the Holland and Powell (2011) model can be applied for the species (in K)
    real Tmax = {};
};

/// The SUPCRTBL parameters for calculating standard thermodynamic properties of minerals without phase transition using Holland and Powell (2011) model.
struct SupcrtblParamsMineralHollandPowell
{
    /// The apparent standard molal Gibbs free energy of formation of the substance from its elements (in kJ/mol)
    real Gf;

    /// The apparent standard molal enthalpy of formation of the substance from its elements (in kJ/mol)
    real Hf;

    /// The standard molal entropy of the substance at reference temperature and pressure (in J/(mol·K))
    real Sr;

    /// The standard molal volume of the substance at reference temperature and pressure (in J/(bar·mol))
    real Vr;

    /// The coefficient `a` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K))
    real a;

    /// The coefficient `b` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K²)) (multiplied by 1e+5)
    real b;

    /// The coefficient `c` of the Holland and Powell (2011) thermodynamic model (in (kJ·K)/mol)
    real c;

    /// The coefficient `d` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K½))
    real d;

    /// The coefficient `α` of the Holland and Powell (2011) thermodynamic model (in 1/K) (multiplied by 1e+5)
    real alpha;

    /// The coefficient `κ` of the Holland and Powell (2011) thermodynamic model (in kbar)
    real kappa;

    /// The coefficient `κ'` of the Holland and Powell (2011) thermodynamic model (dimensionless)
    real kappap;

    /// The coefficient `κ''` of the Holland and Powell (2011) thermodynamic model (in 1/kbar)
    real kappapp;

    /// The number of atoms in the chemical formula of the mineral
    real numatoms;

    /// The maximum temperature at which the Holland and Powell (2011) model can be applied for the species (in K)
    real Tmax = {};
};

/// The SUPCRTBL parameters for calculating standard thermodynamic properties of minerals with phase transition using Holland and Powell (2011) model and Landau theory.
struct SupcrtblParamsMineralHollandPowellLandau
{
    /// The apparent standard molal Gibbs free energy of formation of the substance from its elements (in kJ/mol)
    real Gf;

    /// The apparent standard molal enthalpy of formation of the substance from its elements (in kJ/mol)
    real Hf;

    /// The standard molal entropy of the substance at reference temperature and pressure (in J/(mol·K))
    real Sr;

    /// The standard molal volume of the substance at reference temperature and pressure (in J/(bar·mol))
    real Vr;

    /// The coefficient `a` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K))
    real a;

    /// The coefficient `b` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K²)) (multiplied by 1e+5)
    real b;

    /// The coefficient `c` of the Holland and Powell (2011) thermodynamic model (in (kJ·K)/mol)
    real c;

    /// The coefficient `d` of the Holland and Powell (2011) thermodynamic model (in kJ/(mol·K½))
    real d;

    /// The coefficient `α` of the Holland and Powell (2011) thermodynamic model (in 1/K) (multiplied by 1e+5)
    real alpha;

    /// The coefficient `κ` of the Holland and Powell (2011) thermodynamic model (in kbar)
    real kappa;

    /// The coefficient `κ'` of the Holland and Powell (2011) thermodynamic model (dimensionless)
    real kappap;

    /// The coefficient `κ''` of the Holland and Powell (2011) thermodynamic model (in 1/kbar)
    real kappapp;

    /// The number of atoms in the chemical formula of the mineral
    real numatoms;

    /// The critical temperature of the mineral at 1 bar (in K)
    real Tcr;

    /// The entropy of disordering of the mineral at the critical temperature above (in J/(mol·K))
    real Smax;

    /// The volume of disordering of the mineral at the critical temperature above (in J/(bar·mol))
    real Vmax;

    /// The maximum temperature at which the Holland and Powell (2011) model can be applied for the species (in K)
    real Tmax = {};
};

} // namespace Reaktoro
