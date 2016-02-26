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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Math/LU.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumState;

/// A type that represents a regularized optimization problem.
struct Regularizer
{
    //=============================================================================================
    // Data related to trivial and linearly dependent constraints.
    // Note that these members do not need to be recomputed with matrix A does not change.
    //=============================================================================================
    // The indices of the variables fixed at the lower bound
    Indices itrivial_variables;

    // The indices of the equality constraints whose participating variables are fixed at the lower bound
    Indices itrivial_constraints;

    // The indices of the non-trivial variables
    Indices inontrivial_variables;

    // The indices of the non-trivial constraints
    Indices inontrivial_constraints;

    /// The permutation matrix used to order linearly independent rows.
    PermutationMatrix P_li;

    /// The number of linearly independent rows.
    Index num_li;

    /// The coefficient matrix `A` used in the last regularization.
    Matrix A_last;

    /// The coefficient matrix `A` without trivial constraints and linearly dependent rows.
    /// If new `A` in `update` equals `A_last`, then there is no need to update `A_star`.
    Matrix A_star;

    /// The right-hand side vector `b` without trivial constraints and linearly dependent rows.
    /// If new `A` in `update` equals `A_last`, then there is no need to update `b_star`.
    Matrix b_star;

    //=============================================================================================
    // Data related to regularization related to round-off errors.
    // These members are calculated from `A_star` and `b_star` whenever the weights change order
    // (e.g., the 10th variable has now higher weight than the 5th variable comparing last call.
    //=============================================================================================
    /// The weights computed for the regularization.
    Vector W;

    // The regularizer matrix that is applied to the coefficient matrix `A_star` as `reg(A) = R*A_star`.
    /// If the new set of basic variables are the same as last, then there is no need to update `R`.
    Matrix R;

    // The inverse of the regularizer matrix R.
    /// If the new set of basic variables are the same as last, then there is no need to update `invR`.
    Matrix invR;

    /// The coefficient matrix computed as `A_reg = R * A_star`.
    /// If the new set of basic variables are the same as last, then there is no need to update `A_reg`.
    Matrix A_reg;

    /// The right-hand side vector `b_reg` computed as `b_reg = R * b_star`.
    /// This member must always be updated because `b` changes frequently.
    Matrix b_reg;

    // The indices of basic/independent variables that compose the others.
    Indices ibasic_variables;

    // The indices of basic/independent variables that compose the others in the last call.
    Indices ibasic_variables_last;

    // The full-pivoting LU decomposition of the coefficient matrix `A_star * diag(W)`.
    LU lu;

    // The evaluation of the original objective function.
    ObjectiveResult f;

    /// The regularized OptimumProblem instance.
    OptimumProblem problem;

    /// The regularized OptimumState instance.
    OptimumState state;

    /// The regularized OptimumOptions instance.
    OptimumOptions options;

    /// Update the regularization of the optimum components.
    auto update(const OptimumProblem& problem, const OptimumState& state, const OptimumOptions& options) -> void;

    /// Regularize the constraints by removing trivial constraints that forces variables to be on the bounds.
    auto regularizeTrivialConstraints() -> void;

    /// Regularize the constraints by removing linearly independent constraints.
    auto regularizeLinearlyDependentConstraints() -> void;

    /// Regularize the constraints by transforming them in a cannonical form which helps preventing round-off errors.
    auto regularizeNonCannonicalConstraints() -> void;

    /// Regularize the constraints by adjusting the constraint right-hand side values if they cause infeasible solutions.
    auto regularizeInfeasibleConstraints() -> void;

    /// Regularize an OptimumState instance.
    /// This method should only be called after the optimization problem has been regularized.
    /// @param orig The original OptimumState instance.
    /// @param[in,out] reg The regularized OptimumState instance.
    auto regularizeOptimumState(const OptimumState& orig, OptimumState& reg) -> void;

    /// Regularize an OptimumOptions instance.
    /// This method should only be called after the optimization problem has been regularized.
    /// @param orig The original OptimumOptions instance.
    /// @param[in,out] reg The regularized OptimumOptions instance.
    auto regularizeOptimumOptions(const OptimumOptions& orig, OptimumOptions& reg) -> void;
};

} // namespace Reaktoro
