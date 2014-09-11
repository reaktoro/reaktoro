/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <functional>
#include <tuple>

// Reaktor includes
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

/**
 * Defines a type that describes the required parameters for aqueous activity calculations
 * @see AqueousActivityParams
 */
struct AqueousActivityParams
{
    /// The temperature for the activity calculation (in units of K)
    double T;

    /// The pressure for the activity calculation (in units of Pa)
    double P;

    /// The molar composition of the aqueous mixture (in units of mol)
    Vector n;

    /// The effective ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ThermoScalar Ie;

    /// The stoichiometric ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ThermoScalar Is;

    /// The molar fractions x of the aqueous species and its molar derivatives
    ThermoVector x;

    /// The molalities of the aqueous species and its molar derivatives (in units of mol/kg)
    ThermoVector m;

    /// The stoichiometric molalities of the ionic species and its molar derivatives (in units of mol/kg)
    ThermoVector ms;

    /// Checks for equality of the aqueous activity parameters
    auto operator==(const AqueousActivityParams& params) const -> bool
    {
        return T == params.T and P == params.P and arma::all(n == params.n);
    }
};

/**
 * Defines the function signature of an aqueous activity function
 * @param params An instance of \ref AqueousActivityParams containing the necessary parameters for the activity calculation
 * @return An instance of @ref ThermoScalar containing the calculated activity and its molar derivatives
 * @see AqueousActivityParams, ThermoScalar
 */
using AqueousActivity = std::function<ThermoScalar(const AqueousActivityParams& params)>;

} // namespace Reaktor
