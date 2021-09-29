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

// Forward declarations
struct StandardThermoProps;

/// The complete set of standard thermodynamic properties of a chemical species.
struct SpeciesThermoProps
{
    /// The temperature used to compute the standard thermodynamic properties (in K).
    real T;

    /// The pressure used to compute the standard thermodynamic properties (in Pa).
    real P;

    /// The standard molar Gibbs energy @f$G^{\circ}@f$ of the species (in J/mol).
    real G0;

    /// The standard molar enthalpy @f$H^{\circ}@f$ of the species (in J/mol).
    real H0;

    /// The standard molar volume @f$V^{\circ}@f$ of the species (in m3/mol).
    real V0;

    /// The temperature derivative of the standard molar volume @f$\partial V^{\circ}/\partial T@f$ of the species (in m³/(mol·K)).
    real VT0;

    /// The pressure derivative of the standard molar volume @f$\partial V^{\circ}/\partial P@f$ of the species (in m³/(mol·Pa)).
    real VP0;

    /// The standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in J/(mol·K)).
    real Cp0;

    /// The standard molar isochoric heat capacity @f$C_{V}^{\circ}@f$ of the species (in J/(mol·K)).
    real Cv0;

    /// The standard molar internal energy @f$U^{\circ}@f$ of the species (in J/mol).
    real U0;

    /// The standard molar entropy @f$S^{\circ}@f$ of the species (in J/(mol·K)).
    real S0;

    /// The standard molar Helmholtz energy @f$A^{\circ}@f$ of the species (in J/mol).
    real A0;

    /// Construct a SpeciesThermoProps object.
    /// @param T The temperature corresponding to the standard thermodynamic properties (in K).
    /// @param P The pressure corresponding to the standard thermodynamic properties (in Pa).
    /// @param sprops The primary standard thermodynamic properties of a chemical species.
    SpeciesThermoProps(const real& T, const real& P, const StandardThermoProps& sprops);
};

} // namespace Reaktoro
