// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Optimization/NonlinearSolver.hpp>

namespace Reaktoro {

/// The options for the description of the Hessian of the Gibbs energy function
enum class GibbsHessian
{
    /// The Hessian of the Gibbs energy function is `H = H(exact)`.
    Exact,

    /// The Hessian of the Gibbs energy function is `H = diag(H(exact))`.
    ExactDiagonal,

    /// The Hessian of the Gibbs energy function is `H = d(ln(x))/dn`, where `x` is the mole fractions of the species.
    Approximation,

    /// The Hessian of the Gibbs energy function is `H = diag(d(ln(x))/dn)`, where `x` is the mole fractions of the species.
    ApproximationDiagonal,
};

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// Construct a default EquilibriumOptions instance
    EquilibriumOptions();

    /// Construct a custom EquilibriumOptions instance by parsing a string
    EquilibriumOptions(const char* str);

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
    GibbsHessian hessian = GibbsHessian::ApproximationDiagonal;

    /// The optimisation method to be used for the equilibrium calculation.
    OptimumMethod method = OptimumMethod::IpNewton;

    /// The options for the optimisation calculation.
    OptimumOptions optimum;

    /// The options for the nonlinear solver.
    NonlinearOptions nonlinear;
};

} // namespace Reaktoro
