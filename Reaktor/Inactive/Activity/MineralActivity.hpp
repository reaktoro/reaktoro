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
#include <Reaktor/Common/PartialScalar.hpp>
#include <Reaktor/Common/PartialVector.hpp>

namespace Reaktor {

/**
 * Defines a type that describes the required parameters for mineral activity calculations
 * @see MineralActivity
 */
struct MineralActivityParams
{
    /// The temperature for the activity calculation (in units of K)
    double T;

    /// The pressure for the activity calculation (in units of Pa)
    double P;

    /// The molar composition of the mineral mixture (in units of mol)
    Vector n;

    /// The molar fractions \b x of all mineral species and their molar derivatives
    PartialVector x;

    /// Checks for equality of the mineral activity parameters
    auto operator==(const MineralActivityParams& params) const -> bool
    {
        return T == params.T and P == params.P and n == params.n;
    }
};

/**
 * Defines the function signature of a mineral activity function
 * @param params An instance of \ref MineralActivityParams containing the necessary parameters for the activity calculation
 * @return An instance of @ref PartialScalar containing the calculated activity and its molar derivatives
 * @see MineralActivityParams, PartialScalar
 */
using MineralActivity = std::function<PartialScalar(const MineralActivityParams& params)>;

} // namespace Reaktor
