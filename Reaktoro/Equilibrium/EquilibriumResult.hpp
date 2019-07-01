// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>

namespace Reaktoro {

/// Provide timing information of the operations during an equilibrium calculation.
struct EquilibriumTiming
{
    /// The time spent for solving the chemical equilibrium problem.
    double solve = 0.0;

    /// The time spent for computing the standard thermochemical properties of the system.
    double standard_thermodynamic_properties = 0.0;

    /// The time spent for computing the chemical properties of the system.
    double chemical_properties = 0.0;

    /// Self addition of another EquilibriumTiming instance to this one.
    auto operator+=(const EquilibriumTiming& other) -> EquilibriumTiming&;
};

/// A type used to describe the result of an equilibrium calculation
struct EquilibriumResult
{
    /// The result of the optimisation calculation.
    OptimumResult optimum;

    /// The timing information of the operations during an equilibrium calculation.
    EquilibriumTiming timing;

    /// Self addition assignment to accumulate results.
    auto operator+=(const EquilibriumResult& other) -> EquilibriumResult&;
};

} // namespace Reaktoro
