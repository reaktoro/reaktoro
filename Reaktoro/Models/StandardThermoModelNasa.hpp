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
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

/// The parameters in the NASA polynomial model for calculating standard thermodynamic properties of gases and condensed species.
struct StandardThermoModelParamsNasa
{
    /// Used to store Nasa polynomial coefficients valid for a certain temperature interval
    struct Polynomial
    {
        double Tmin; ///< The minimum temperature (in K) for which the NASA polynomial coefficients below are valid.
        double Tmax; ///< The maximum temperature (in K) for which the NASA polynomial coefficients below are valid.

        String label; ///< The label describing the species name and its physical form or crystal configuration such as Mg(L)L, Li3ALF6(IV), etc. (optional).
        AggregateState state; ///< The aggregate state of the species within this temperature interval.

        Param a1; /// The least-square coefficient \eq{a_1} in the regression model for \eq{C_{p}^{\circ}}.
        Param a2; /// The least-square coefficient \eq{a_2} in the regression model for \eq{C_{p}^{\circ}}.
        Param a3; /// The least-square coefficient \eq{a_3} in the regression model for \eq{C_{p}^{\circ}}.
        Param a4; /// The least-square coefficient \eq{a_4} in the regression model for \eq{C_{p}^{\circ}}.
        Param a5; /// The least-square coefficient \eq{a_5} in the regression model for \eq{C_{p}^{\circ}}.
        Param a6; /// The least-square coefficient \eq{a_6} in the regression model for \eq{C_{p}^{\circ}}.
        Param a7; /// The least-square coefficient \eq{a_7} in the regression model for \eq{C_{p}^{\circ}}.
        Param b1; /// The integration constant \eq{b_1} used to compute \eq{H^\circ}.
        Param b2; /// The integration constant \eq{b_2} used to compute \eq{S^\circ}.
    };

    /// The heat of formation \eq{\Delta H_{f}^{\circ}} at 298.15 K (in J/mol).
    real dHf;

    /// The value of \eq{\Delta H_{0}^{\circ}=H^{\circ}(298.15)-H^{\circ}(0)} (in J/mol).
    real dH0;

    /// The assigned enthalpy (in J/mol) of the species when there are no temperature intervals.
    real H0;

    /// The temperature (in K) corresponding to the assigned enthalpy when there are no temperature intervals.
    double T0;

    /// The Nasa polynomials for one or more temperature intervals.
    Vec<Polynomial> polynomials;
};

/// Return a function that calculates thermodynamic properties of a species using the Maier-Kelley model.
auto StandardThermoModelNasa(const StandardThermoModelParamsNasa& params) -> StandardThermoModel;

//=================================================================================================
// AUXILIARY METHODS
//=================================================================================================

namespace detail {

/// Return the index of the temperature interval that contains a given temperature.
/// If the given temperature is out-of-bounds, then the number of intervals is returned.
/// @param polynomials The NASA polynomials for one or more temperature intervals
/// @param T The temperature (in K)
auto indexTemperatureInterval(const Vec<StandardThermoModelParamsNasa::Polynomial>& polynomials, const real& T) -> Index;

/// Compute the standard thermodynamic properties of a species with given NASA thermodynamic parameters.
/// @param polynomial The NASA polynomial covering a certain temperature interval
/// @param T The temperature (in K)
auto computeStandardThermoProps(const StandardThermoModelParamsNasa::Polynomial& polynomial, const real& T) -> StandardThermoProps;

/// Compute the standard thermodynamic properties of a species with given NASA thermodynamic parameters.
/// @param params The parameters in the NASA polynomial model for a species
/// @param T The temperature (in K)
auto computeStandardThermoProps(const StandardThermoModelParamsNasa& params, const real& T) -> StandardThermoProps;

} // namespace detail
} // namespace Reaktoro
