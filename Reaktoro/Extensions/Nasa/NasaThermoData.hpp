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

/// The thermodynamic data for a substance in a NASA database.
/// The coefficients and integration constants listed here are documented in
/// [1] (Section 4.2, page 19) and also in [2] (Appendix A, page 74).
///
/// 1. Mcbride, B.J., Gordon, S., Mcbride, B.J. (1994). Computer program for
///    calculation of complex chemical equilibrium compositions and applications.
///    I: Analysis. In NASA Reference Publication 1311. https://doi.org/NASA RP-1311
/// 2. McBride, B.J., Gordon, S. (1996). Computer Program for Calculation of
///    Complex Chemical Equilibrium Compositions and Applications: II-User Manual
///    and Program Description. In NASA Reference Publication 1311.
///    https://doi.org/NASA RP-1311
struct NasaThermoData
{
    double Tmin = {}; ///< The minimum temperature (in K) for which the coefficients below are valid.
    double Tmax = {}; ///< The maximum temperature (in K) for which the coefficients below are valid.

    long qN = 7; ///< The number of exponent coefficients \eq{q_i}, which is always 7.

    double q1 = -2.0; ///< The exponent \eq{q_1} in the regression model for \eq{C_{p}^{\circ}}.
    double q2 = -1.0; ///< The exponent \eq{q_2} in the regression model for \eq{C_{p}^{\circ}}.
    double q3 =  0.0; ///< The exponent \eq{q_3} in the regression model for \eq{C_{p}^{\circ}}.
    double q4 =  1.0; ///< The exponent \eq{q_4} in the regression model for \eq{C_{p}^{\circ}}.
    double q5 =  2.0; ///< The exponent \eq{q_5} in the regression model for \eq{C_{p}^{\circ}}.
    double q6 =  3.0; ///< The exponent \eq{q_6} in the regression model for \eq{C_{p}^{\circ}}.
    double q7 =  4.0; ///< The exponent \eq{q_7} in the regression model for \eq{C_{p}^{\circ}}.

    real a1; ///< The least-square coefficient \eq{a_1} in the regression model for \eq{C_{p}^{\circ}}.
    real a2; ///< The least-square coefficient \eq{a_2} in the regression model for \eq{C_{p}^{\circ}}.
    real a3; ///< The least-square coefficient \eq{a_3} in the regression model for \eq{C_{p}^{\circ}}.
    real a4; ///< The least-square coefficient \eq{a_4} in the regression model for \eq{C_{p}^{\circ}}.
    real a5; ///< The least-square coefficient \eq{a_5} in the regression model for \eq{C_{p}^{\circ}}.
    real a6; ///< The least-square coefficient \eq{a_6} in the regression model for \eq{C_{p}^{\circ}}.
    real a7; ///< The least-square coefficient \eq{a_7} in the regression model for \eq{C_{p}^{\circ}}.
    real b1; ///< The integration constant \eq{b_1} used to compute \eq{H^\circ}.
    real b2; ///< The integration constant \eq{b_2} used to compute \eq{S^\circ}.
};

/// Return true if two NasaThermoData objects are different.
auto operator!=(const NasaThermoData& l, const NasaThermoData& r) -> bool;

/// Return true if two NasaThermoData objects are equal.
auto operator==(const NasaThermoData& l, const NasaThermoData& r) -> bool;

} // namespace Reaktoro
