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
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

/// The thermodynamic properties of water computed in PHREEQC.
struct PhreeqcWaterThermoProps
{
    /// The density of water (in g/cm3).
    real rho_0 = {};

    /// The compressibility of water, d(ln(rho))/dP, (in 1/atm).
    real kappa_0 = {};
};

/// The electrostatic properties of water computed in PHREEQC.
struct PhreeqcWaterElectroProps
{
    /// The relative dielectric constant of water.
    real eps_r = {};

    // The Debye-Hueckel A coefficient (in (mol/kg)^-0.5).
    real DH_A = {};

    // The Debye-Hueckel B coefficient (in 1/Angstrom(mol/kg)^-0.5).
    real DH_B = {};

    // The Debye-Hueckel limiting slope (in (cm3/mol)(mol/kg)^-0.5).
    real DH_Av = {};

    // The Born function (-1/eps_r + 1) * 41.84004, for supcrt calc'n of molal volume
    real ZBrn = {};

    // The Born function d(ln(eps_r))/dP / eps_r * 41.84004, for supcrt calc'n of molal volume
    real QBrn = {};
};

/// The thermodynamic and electrostatic properties of water computed in PHREEQC.
struct PhreeqcWaterProps
{
    /// The thermodynamic properties of water.
    PhreeqcWaterThermoProps wtp;

    /// The electrostatic properties of water.
    PhreeqcWaterElectroProps wep;
};

/// Compute the thermodynamic properties of water using same model as used in PHREEQC.
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
auto waterThermoProps(real T, real P) -> PhreeqcWaterThermoProps;

/// Compute the electrostatic properties of water using same model as used in PHREEQC.
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
/// @param wtp The thermodynamic properties of water at given temperature and pressure
auto waterElectroProps(real T, real P, PhreeqcWaterThermoProps wtp) -> PhreeqcWaterElectroProps;

/// Compute the thermodynamic and electrostatic properties of water using same model as used in PHREEQC.
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
auto waterProps(real T, real P) -> PhreeqcWaterProps;

} // namespace PhreeqcUtils
} // namespace Reaktoro
