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

// C++ includes
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

/// The complete set of standard thermodynamic properties of a chemical reaction.
struct ReactionThermoProps
{
    /// The temperature used to compute the reaction properties (in K).
    real T;

    /// The pressure used to compute the reaction properties (in Pa).
    real P;

    /// The equilibrium constant of the reaction (log base 10).
    real lgK;

    /// The change in standard molar Gibbs energy @f$\Delta G^{\circ}@f$ in the reaction (in J/mol).
    real dG0;

    /// The change in standard molar enthalpy @f$\Delta H^{\circ}@f$ in the reaction (in J/mol).
    real dH0;

    /// The change in standard molar volume @f$\Delta V^{\circ}@f$ in the reaction (in m3/mol).
    real dV0;

    /// The change in temperature derivative of the standard molar volume @f$\partial \Delta V^{\circ}/\partial T@f$ in the reaction (in m³/(mol·K)).
    real dVT0;

    /// The change in pressure derivative of the standard molar volume @f$\partial \Delta V^{\circ}/\partial P@f$ in the reaction (in m³/(mol·Pa)).
    real dVP0;

    /// The change in standard molar isobaric heat capacity @f$\Delta C_{P}^{\circ}@f$ in the reaction (in J/(mol·K)).
    real dCp0;

    /// The change in standard molar isochoric heat capacity @f$\Delta C_{V}^{\circ}@f$ in the reaction (in J/(mol·K)).
    real dCv0;

    /// The change in standard molar internal energy @f$\Delta U^{\circ}@f$ in the reaction (in J/mol).
    real dU0;

    /// The change in standard molar entropy @f$\Delta S^{\circ}@f$ in the reaction (in J/(mol·K)).
    real dS0;

    /// The change in standard molar Helmholtz energy @f$\Delta A^{\circ}@f$ in the reaction (in J/mol).
    real dA0;
};

/// Output a ReactionThermoProps object to an output stream.
auto operator<<(std::ostream& out, const ReactionThermoProps& props) -> std::ostream&;

} // namespace Reaktoro
