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
#include <Reaktoro/Core/ReactionThermoProps.hpp>

namespace Reaktoro {

/// Return a function that calculates thermodynamic properties of a reaction using van't Hoff's model.
///
/// In this model, the equilibrium constant of the reaction
/// is computed at a given temperature with:
///
/// @eqc{\ln K=\ln K_{0}-\frac{\Delta H_{0}^{\circ}}{R}\left(\frac{1}{T}-\frac{1}{T_{0}}\right),}
///
/// where @eq{\ln K_{0}} and @eq{\Delta H_{0}^{\circ}} are the equilibrium constant and
/// enthalpy of reaction at reference temperature @eq{T_0}. From this model, we can calculate
/// the standard Gibbs energy of reaction using:
///
/// @eqc{\Delta G^{\circ}=-RT\ln K}
///
/// while the standard enthalpy of reaction is assumed constant:
///
/// @eqc{\Delta H^{\circ}=\Delta H_{0}^{\circ}.}
///
/// @param lgK0 The equilibrium constant (log base 10) of the reaction at reference temperature.
/// @param dH0 The change in standard molar enthalpy of the reaction at reference temperature (in J/mol).
/// @param Tr The reference temperature (in K).
auto ReactionThermoModelVantHoff(Param lgK0, Param dH0, Param Tr) -> ReactionThermoModel;

} // namespace Reaktoro
