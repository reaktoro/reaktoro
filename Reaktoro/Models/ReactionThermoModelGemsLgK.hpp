// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

/// The parameters in the thermodynamic model for a formation reaction based on GEM-Selektor's expression for @eq{\lg K(T)}.
struct ReactionThermoModelParamsGemsLgK
{
    /// The coefficient @eq{A_0} in the reaction thermodynamic model.
    Param A0;

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

/// Return a function that calculates thermodynamic properties of a reaction using GEM-Selektor's expression for @eq{\lg K(T)}.
///
/// In this model, the equilibrium constant of the reaction is computed at a given temperature with:
///
/// @eqc{\log_{10}K=A_{0}+A_{1}T+A_{2}T^{-1}+A_{3}\ln T+A_{4}T^{-2}+A_{5}T^{2}+A_{6}T^{-0.5},}
///
/// where @eq{A_i} are given coefficients. From this model, we can calculate
/// the standard Gibbs energy of reaction using:
///
/// @eqc{\Delta G^{\circ}=-RT\left(A_{0}+A_{1}T+A_{2}T^{-1}+A_{3}\ln T+A_{4}T^{-2}+A_{5}T^{2}+A_{6}T^{-0.5}\right)\ln_{10}}
///
/// and the standard enthalpy of reaction using:
///
/// @eqc{\Delta H^{\circ}=R\left(A_{1}T^{2}-A_{2}+A_{3}T-2A_{4}T^{-1}+2A_{5}T^{3}-0.5A_{6}T^{0.5}\right)\ln_{10},}
///
/// considering that:
///
/// @eqc{\Delta G^{\circ}=-RT\ln K}
///
/// and
///
/// @eqc{\Delta H^{\circ}\equiv T^{2}\frac{\partial}{\partial T}\left(-\frac{\Delta G^{\circ}}{T}\right)=RT^{2}\frac{\partial\ln K}{\partial T}.}
///
/// Reference:
/// - Kulik, D., Wagner T. (2010) Part 3. Temperature corrections of
///   standard molar (partial molal) thermodynamic properties of substances and
///   reactions using data in ReacDC records of GEM-Selektor.
/// - Haas, J.L., Fisher, J.R. (1976). Simultaneous evaluation and correlation
///   of thermodynamic data. In American Journal of Science (Vol. 276, Issue 4,
///   pp. 525–545). https://doi.org/10.2475/ajs.276.4.525
auto ReactionThermoModelGemsLgK(const ReactionThermoModelParamsGemsLgK& params) -> ReactionThermoModel;

} // namespace Reaktoro
