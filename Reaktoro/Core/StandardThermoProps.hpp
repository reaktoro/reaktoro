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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Model.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic properties of a chemical species.
struct StandardThermoProps
{
    /// The standard molar Gibbs energy @f$G^{\circ}@f$ of formation of the species (in J/mol).
    real G0;

    /// The standard molar enthalpy @f$H^{\circ}@f$ of formation of the species (in J/mol).
    real H0;

    /// The standard molar volume @f$V^{\circ}@f$ of the species (in m³/mol).
    real V0;

    /// The standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in J/(mol·K)).
    real Cp0;

    /// The temperature derivative of the standard molar volume @f$\partial V^{\circ}/\partial T@f$ of the species (in m³/(mol·K)).
    real VT0;

    /// The pressure derivative of the standard molar volume @f$\partial V^{\circ}/\partial P@f$ of the species (in m³/(mol·Pa)).
    real VP0;
};

/// The function type for calculation of standard thermodynamic properties of a species.
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
/// @return The standard thermodynamic properties of the species
using StandardThermoModel = Model<StandardThermoProps(real T, real P)>;

} // namespace Reaktoro
