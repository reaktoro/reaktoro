// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/ReactionStandardThermoModel.hpp>

namespace Reaktoro {

/// The parameters in a thermodynamic model for a formation reaction based on constant @eq{\lg K(T)}.
struct ReactionStandardThermoModelParamsConstLgK
{
    /// The equilibrium constant @eq{\lg K_{\mathrm{r}}} (log base 10) of the reaction at @eq{T_{\mathrm{r}}} and @eq{P_{\mathrm{r}}}.
    real lgKr;

    /// The reference pressure @eq{P_{\mathrm{r}}} (in Pa).
    real Pr = 100'000;
};

/// Return a function that calculates thermodynamic properties of a reaction using a constant model for @eq{\lg K(T)}.
///
/// In this model, the equilibrium constant of the reaction remains constant:
///
/// @eqc{\lg K(T)=\lg K_{\mathrm{r}},}
///
/// where @eq{\lg K_{\mathrm{r}}} is a given equilibrium constant at a reference temperature and pressure.
///
/// The standard Gibbs energy of the reaction is computed using:
///
/// @eqc{\Delta G^{\circ}(T,P)=-RT\ln K_{\mathrm{r}}+\Delta V^{\circ}(P-P_{\mathrm{r}}),}
///
/// the standard enthalpy of the reaction using:
///
/// @eqc{\Delta H^{\circ}(T,P)=\Delta V^{\circ}(P-P_{\mathrm{r}}).}
///
/// and the standard isobaric heat capacity of reaction using:
///
/// @eqc{\Delta C_{P}^{\circ}(T, P) = \Delta V^{\circ}_T(P-P_{\mathrm{r}}),}
///
/// where we have considered the following thermodynamic relations:
///
/// @eqc{\begin{align*} \Delta G^{\circ} & \equiv-RT\ln K,\vphantom{\left(-\frac{\Delta G^{\circ}}{T}\right)}\\ \Delta H^{\circ} & \equiv T^{2}\frac{\partial}{\partial T}\left(-\frac{\Delta G^{\circ}}{T}\right)_{P}=RT^{2}\frac{\partial\ln K}{\partial T},\\ \Delta C_{P}^{\circ} & \equiv\left(\dfrac{\partial\Delta H^{\circ}}{\partial T}\right)_{P}. \end{align*}}
///
/// Note that a pressure correction is introduced above, where @eq{\Delta V^{\circ}}
/// is the change of standard molar volume of the reaction at @eq{T} and @eq{P},
/// with @eq{P_{\mathrm{r}}} denoting a given reference pressure. @eq{\Delta V^{\circ}_T}
/// is the temperature derivative of @eq{\Delta V^{\circ}} at constant pressure.
auto ReactionStandardThermoModelConstLgK(const ReactionStandardThermoModelParamsConstLgK& params) -> ReactionStandardThermoModel;

} // namespace Reaktoro
