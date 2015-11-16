// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumMethod.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>

namespace Reaktoro {

/// The options for the description of the Hessian of the Gibbs energy function
enum class GibbsHessian
{
    /// The Hessian approximation given by the diagonal matrix `H = diag(inv(n))`.
    Diagonal,

    /// The Hessian approximation given by the diagonal matrix `H(i,i) = 1/ni` if species
    /// `i` lives in a multicomponent phase, otherwise, `H(i,i) = 0`.
    SparseDiagonal,

    /// The exact Hessian of the Gibbs energy function.
    Exact
};

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// The parameter ε for the numerical representation of a zero molar amount.
    /// The molar amount of the `i`-th species is considered zero if `n[i] < ε*min(b)`,
    /// where `b` is the vector of element molar amounts.
    double epsilon = 1e-20;

    /// The boolean flag that indicates if warm-start strategy should be used when possible.
    /// Setting this flag to true will cause equilibrium calculations to use the currect
    /// chemical state as an initial guess to the equilibrium calculation. If the current
    /// chemical state is detected to be uninitiallized (e.g., all species with
    /// zero molar amounts), then cold-start is inevitably used. In this case, a first estimate
    /// of a initial guess will be done using a simplex algorithm, which in most cases generates
    /// a chemical state that works well as initial guess for all equilibrium algorithms.
    bool warmstart = true;

    /// The calculation mode of the Hessian of the Gibbs energy function
    GibbsHessian hessian = GibbsHessian::SparseDiagonal;

    /// The optimisation method to be used for the equilibrium calculation.
    OptimumMethod method = OptimumMethod::IpNewton;

    /// The options for the optimisation calculation.
    OptimumOptions optimum;
};

} // namespace Reaktoro
