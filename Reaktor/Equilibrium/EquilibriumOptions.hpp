// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>

namespace Reaktor {

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// The parameter ε for the numerical representation of a zero molar amount.
    /// The molar amount of the `i`-th species is considered zero if `n[i] < ε*min(b)`,
    /// where `b` is the vector of element molar amounts.
    double epsilon = 1e-50;

    /// The factor τ for the barrier parameter μ defined here as μ = ετ.
    /// The parameter ε is the numerical zero for a molar amount.
    /// @see epsilon
    double tau = 1e-5;

    /// The options for the optimisation calculation.
    OptimumOptions optimum;
};

} // namespace Reaktor
