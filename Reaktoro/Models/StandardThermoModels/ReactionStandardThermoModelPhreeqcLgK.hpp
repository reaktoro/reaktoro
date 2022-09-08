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
#include <Reaktoro/Core/ReactionStandardThermoModel.hpp>

namespace Reaktoro {

/// The parameters in the thermodynamic model for a formation reaction based on PHREEQC's expression for @eq{\lg K(T)}.
struct ReactionStandardThermoModelParamsPhreeqcLgK
{
    /// The coefficient @eq{A_1} in the reaction thermodynamic model.
    Param A1;

    /// The coefficient @eq{A_2} in the reaction thermodynamic model.
    Param A2;

    /// The coefficient @eq{A_3} in the reaction thermodynamic model.
    Param A3;

    /// The coefficient @eq{A_4} in the reaction thermodynamic model.
    Param A4;

    /// The coefficient @eq{A_5} in the reaction thermodynamic model.
    Param A5;

    /// The coefficient @eq{A_6} in the reaction thermodynamic model.
    Param A6;

    /// The reference pressure @eq{P_{\mathrm{r}}} (in Pa).
    real Pr = 100'000;
};

/// Return a function that calculates thermodynamic properties of a reaction using PHREEQC's analytical expression.
///
/// In this model, the equilibrium constant of the reaction is computed at a given temperature with:
///
/// @eqc{\lg K(T)=A_{1}+A_{2}T+A_{3}T^{-1}+A_{4}\log_{10}T+A_{5}T^{-2}+A_{6}T^{2},}
///
/// where @eq{A_i} are given coefficients.
///
/// The standard Gibbs energy of reaction is computed using:
///
/// @eqc{\Delta G^{\circ}(T,P)=-RT\left(A_{1}+A_{2}T+A_{3}T^{-1}+A_{4}\log_{10}T+A_{5}T^{-2}+A_{6}T^{2}\right)\ln_{10}+\Delta V^{\circ}(P-P_{\mathrm{r}}),}
///
/// the standard enthalpy of reaction using:
///
/// @eqc{\Delta H^{\circ}(T,P)=R\left(A_{2}T^{2}-A_{3}+\frac{A_{4}}{\ln10}T-2A_{5}T^{-1}+2A_{6}T^{3}\right)\ln_{10}+\Delta V^{\circ}(P-P_{\mathrm{r}}),}
///
/// and the standard isobaric heat capacity of reaction using:
///
/// @eqc{\Delta C_{P}^{\circ}(T,P)=R\left(2A_{2}T+\frac{A_{4}}{\ln10}+2A_{5}T^{-2}+6A_{6}T^{2}\right)\ln_{10}+\Delta V_{T}^{\circ}(P-P_{\mathrm{r}}),}
///
/// where we have considered the following thermodynamic relations:
///
/// @eqc{\begin{align*} \Delta G^{\circ} & \equiv-RT\ln K,\vphantom{\left(-\frac{\Delta G^{\circ}}{T}\right)}\\ \Delta H^{\circ} & \equiv T^{2}\frac{\partial}{\partial T}\left(-\frac{\Delta G^{\circ}}{T}\right)_{P}=RT^{2}\frac{\partial\ln K}{\partial T},\\ \Delta C_{P}^{\circ} & \equiv\left(\dfrac{\partial\Delta H^{\circ}}{\partial T}\right)_{P}. \end{align*}}
///
/// Note that a pressure correction is introduced above, where @eq{\Delta V^{\circ}}
/// is the change of standard molar volume of the reaction at @eq{T} and @eq{P},
/// with @eq{P_{\mathrm{r}}} denoting a given reference pressure. @eq{\Delta V^{\circ}_T}
/// is the temperature derivative of @eq{\Delta V^{\circ}} at constant pressure.
///
/// Reference:
/// - Parkhurst, D.L., Appelo, C.A.J. (2013). Description of input and examples
///   for PHREEQC version 3—A computer program for speciation, batch-reaction,
///   one-dimensional transport, and inverse geochemical calculations. In
///   Groundwater Book 6, Modeling Techniques (p. 497). U.S. Geological Survey
///   Techniques and Methods. http://pubs.usgs.gov/tm/06/a43
auto ReactionStandardThermoModelPhreeqcLgK(const ReactionStandardThermoModelParamsPhreeqcLgK& params) -> ReactionStandardThermoModel;

} // namespace Reaktoro
