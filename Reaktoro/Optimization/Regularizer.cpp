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

#include "Regularizer.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Math/LU.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>

namespace Reaktoro {

struct Regularizer::Impl
{
	/// The parameters for the regularization
	OptimumParamsRegularization params;

	/// The evaluation of the original objective function.
	ObjectiveResult f;

	//=============================================================================================
	// Data related to trivial and linearly dependent constraints.
	// Note that these members do not need to be recomputed with matrix A does not change.
	//=============================================================================================
	/// The indices of the variables fixed at the lower bound
	Indices itrivial_variables;

	/// The indices of the equality constraints whose participating variables are fixed at the lower bound
	Indices itrivial_constraints;

	/// The indices of the non-trivial variables
	Indices inontrivial_variables;

	/// The indices of the non-trivial constraints
	Indices inontrivial_constraints;

	/// The indices of the linearly independent constraints
	Indices ili_constraints;

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
	// Data related to echelonization of the constraints (helps with round-off errors).
	// These members are calculated from `A_star` and `b_star`.
	//=============================================================================================
	/// The weights computed for the regularization.
	Vector W;

	/// The regularizer matrix that is applied to the coefficient matrix `A_star` as `reg(A) = R*A_star`.
	/// If the new set of basic variables are the same as last, then there is no need to update `R`.
	Matrix R;

	/// The inverse of the regularizer matrix R.
	/// If the new set of basic variables are the same as last, then there is no need to update `invR`.
	Matrix invR;

	/// The coefficient matrix computed as `A_echelon = R * A_star`.
	/// If the new set of basic variables are the same as last, then there is no need to update `A_echelon`.
	Matrix A_echelon;

	/// The right-hand side vector `b_echelon` computed as `b_echelon = R * b_star`.
	/// This member must always be updated because `b` changes frequently.
	Matrix b_echelon;

	// The indices of basic/independent variables that compose the others.
	Indices ibasic_variables;

	/// The indices of basic/independent variables that compose the others in the last call.
	Indices ibasic_variables_last;

	/// The full-pivoting LU decomposition of the coefficient matrix `A_star * diag(W)`.
	LU lu;

    /// Determine the trivial constraints and trivial variables.
    /// Trivial constraints are all those which fix the values of
    /// some variables (trivial variables) to the bounds.
    auto determineTrivialConstraints(const OptimumProblem& problem) -> void;

    /// Determine the linearly dependent constraints.
    /// This method should be called only after `determineTrivialConstraints`.
    auto determineLinearlyDependentConstraints(const OptimumProblem& problem) -> void;

    /// Assemble the constraints in cannonical form to help in the prevention of round-off errors.
    /// This method should be called only after `determineLinearlyDependentConstraints`.
    auto assembleEchelonConstraints(const OptimumState& state) -> void;

    /// Remove all trivial constraints from the optimum problem.
    /// This method also adjust the members in optimum state and options that
    /// correspond to the removed trivial constraints and trivial variables.
    /// This method should be called after `assembleEchelonConstraints`
    auto removeTrivialConstraints(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Remove all linearly dependent constraints from the optimum problem.
    /// This method also adjust the members in optimum state and options that
    /// correspond to the removed linear constraints.
    /// This method should be called after `removeTrivialConstraints`
    auto removeLinearlyDependentConstraints(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Fix the infeasible linear constraints by adjusting their right-hand side values.
    /// This method should be called after `removeLinearlyDependentConstraints`
    auto fixInfeasibleConstraints(OptimumProblem& problem) -> void;

    /// Transform the original linear constraints into a echelon form as a way to minimize round-off errors.
    /// This method should be called after `fixInfeasibleConstraints`
    auto echelonizeConstraints(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Update the linear constraints of an optimum problem.
    /// This is the last method to be called during the regularization steps.
    /// It should be called after `echelonizeConstraints`.
    auto updateConstraints(OptimumProblem& problem) -> void;

    /// Regularize the optimum problem, state, and options before they are used in an optimization calculation.
    auto regularize(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Recover an optimum state to an state that corresponds to the original optimum problem.
    auto recover(const OptimumProblem& problem, OptimumState& state) -> void;
};

auto Regularizer::Impl::determineTrivialConstraints(const OptimumProblem& problem) -> void
{
	// Skip if the new coefficient matrix is the same as last time.
	if(problem.A == A_last)
		return;

    // Auxiliary references
    const Matrix& A = problem.A;
    const Vector& b = problem.b;
    const Vector& l = problem.l;

    // The number of rows and cols in the original coefficient matrix
    const Index m = A.rows();
    const Index n = A.cols();

	// Clear previous states of trivial and non-trivial constraints and variables
	itrivial_constraints.clear();
	itrivial_variables.clear();
	inontrivial_constraints.clear();
	inontrivial_variables.clear();

	// Return true if the i-th constraint forces the variables to be fixed on the lower bounds
	auto istrivial = [&](Index irow)
	{
		return ( min(A.row(irow)) >= 0 && A.row(irow)*l >= b[irow] ) ||
			   ( max(A.row(irow)) <= 0 && A.row(irow)*l <= b[irow] );
	};

	// Determine the original equality constraints that fix variables on the lower bound
	for(Index i = 0; i < m; ++i)
		if(istrivial(i))
			itrivial_constraints.push_back(i);

	// Skip the rest if there are no trivial constraints
	if(itrivial_constraints.size())
	{
		// Determine the original trivial variables that are fixed at their lower bounds
		for(Index i = 0; i < n; ++i)
			for(Index j = 0; j < m; ++j)
				if(A(j, i) != 0.0 && contained(j, itrivial_constraints))
					{ itrivial_variables.push_back(i); break; }

		// Update  the indices of the non-trivial original constraints
		inontrivial_constraints = difference(range(m), itrivial_constraints);

		// Update the indices of the non-trivial original variables
		inontrivial_variables = difference(range(n), itrivial_variables);

		// Assert there not all contraints are trivial
		Assert(inontrivial_variables.size(),
			"Could not accept the optimization problem.",
			"The provided problem contains only trivial constraints.");

		// Initialize the matrix `A_star` by removing trivial constraints and variables
		A_star = submatrix(A, inontrivial_constraints, inontrivial_variables);
	}
	else
	{
		// Initialize the matrix `A_star` as the original matrix A
		A_star = A;
	}
}

auto Regularizer::Impl::determineLinearlyDependentConstraints(const OptimumProblem& problem) -> void
{
	// Skip if the new coefficient matrix is the same as last time.
	if(problem.A == A_last)
		return;

    // The number of rows and cols in the coefficient matrix A*,
	// i.e., the original A matrix with removed trivial constraints and variables
    const Index m = A_star.rows();
    const Index n = A_star.cols();

    // Compute the LU decomposition of A*
    lu.compute(A_star);

    // Auxiliary references to LU components
    const auto& P = lu.P;
    const auto& rank = lu.rank;

    // Check if there are linearly dependent rows in A*
    if(rank != m)
    {
        // Initialize the indices of the linearly independent constraints
    	ili_constraints = Indices(P.indices().data(), P.indices().data() + rank);

        // Update the permutation matrix and the number of linearly independent constraints
        P_li = lu.P;
        num_li = lu.rank;

        // Permute the rows of A and remove the linearly dependent ones
        A_star = P_li * A_star;
        A_star.conservativeResize(num_li, n);
    }
    else
    {
        // There is no need to update A_star here because all rows are linearly independent
    	ili_constraints = range(m);
        P_li.setIdentity(m);
        num_li = m;
    }
}

auto Regularizer::Impl::assembleEchelonConstraints(const OptimumState& state) -> void
{
    // Initialize the weight vector `W`
    W = abs(state.x);

    // The threshold used to avoid scaling by very tiny components (prevent it from being zero)
    const double threshold = 1e-10 * (max(W) + 1);

    // Skip the echelonization if all weights are less than the threshold
    if(max(W) <= threshold)
    	{ A_echelon.conservativeResize(0, 0); return; }

    // Remove very tiny values in W
    W = (W.array() > threshold).select(W, threshold);

    // Compute the LU decomposition of matrix A_star with column-sorting weights.
    // Columsn corresponding to variables with higher weights are moved to the beginning of the matrix.
    // This gives preference for those variables to become linearly independent basis
    lu.compute(A_star, W);

    // Auxiliary references to LU components
    const auto& Q = lu.Q;
    const auto& rank = lu.rank;

    // Initialize the indices of the basic variables
    ibasic_variables = Indices(Q.indices().data(), Q.indices().data() + rank);

    // Check if the new set of basic variables is diffent than the previous
    if(!contained(ibasic_variables, ibasic_variables_last))
    {
        // The rank of the original coefficient matrix
        const auto r = lu.rank;

        // The L factor of the original coefficient matrix
        const auto L = lu.L.topLeftCorner(r, r).triangularView<Eigen::Lower>();

        // The U1 part of U = [U1 U2]
        const auto U1 = lu.U.topLeftCorner(r, r).triangularView<Eigen::Upper>();

        // Compute the regularizer matrix R = inv(U1)*inv(L)
        R = identity(r, r);
        R = L.solve(R);
        R = U1.solve(R);

        // Compute the inverse of the regularizer matrix inv(R) = L*U1
        invR = U1;
        invR = L * invR;

        // Compute the equality constraint regularization
        A_echelon = R * A_star;

        // Check if the regularizer matrix is composed of rationals.
        // If so, round-off errors can be eliminated
        if(params.max_denominator)
        {
            cleanRationalNumbers(A_echelon, params.max_denominator);
            cleanRationalNumbers(R, params.max_denominator);
            cleanRationalNumbers(invR, params.max_denominator);
        }
    }
}

auto Regularizer::Impl::removeTrivialConstraints(
	OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
{
	// Skip the rest if there is no trivial constraint
	if(itrivial_constraints.empty())
		return;

	// Remove trivial components from problem.b
	problem.b = rows(problem.b, inontrivial_constraints);

	// Remove trivial components from problem.c
	if(problem.c.rows())
		problem.c = rows(problem.c, inontrivial_variables);

	// Remove trivial components from problem.l
	if(problem.l.rows())
		problem.l = rows(problem.l, inontrivial_variables);

	// Remove trivial components from problem.u
	if(problem.u.rows())
		problem.u = rows(problem.u, inontrivial_variables);

	// Remove trivial components from problem.objective
	if(problem.objective)
	{
		// The auxiliary vector used in the lambda functions below. The use of this vector `x`
		// ensures that trivial components remain unchanged on the lower bounds.
		Vector x = problem.l;

		// The objective function before it is regularized.
		ObjectiveFunction original_objective = problem.objective;

		// The evaluation of the regularized objective function
		ObjectiveResult res;

		problem.objective = [=](const Vector& X) mutable
		{
			rows(x, inontrivial_variables) = X;

			f = original_objective(x);

			res.val = f.val;
			res.grad = rows(f.grad, inontrivial_variables);
			res.hessian.mode = f.hessian.mode;

			if(f.hessian.dense.size())
				res.hessian.dense = submatrix(f.hessian.dense, inontrivial_variables, inontrivial_variables);
			if(f.hessian.diagonal.size())
				res.hessian.diagonal = rows(f.hessian.diagonal, inontrivial_variables);
			if(f.hessian.inverse.size())
				res.hessian.inverse = submatrix(f.hessian.inverse, inontrivial_variables, inontrivial_variables);

			return res;
		};
	}

    // Remove trivial components corresponding to trivial variables and trivial constraints
    state.x = rows(state.x, inontrivial_variables);
    state.y = rows(state.y, inontrivial_constraints);
    state.z = rows(state.z, inontrivial_variables);

    // Update the names of the constraints and variables accordingly
    if(options.output.active)
    {
        options.output.xnames = extract(options.output.xnames, inontrivial_variables);
        options.output.ynames = extract(options.output.ynames, inontrivial_constraints);
        options.output.znames = extract(options.output.znames, inontrivial_variables);
    }
}

auto Regularizer::Impl::removeLinearlyDependentConstraints(
	OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
{
	// The number of rows of the A* matrix
	const Index m = A_star.rows();

	// Skip the rest if A* has all rows linearly independent
	if(num_li == m)
		return;

    // Remove the components in b corresponding to linearly dependent constraints
    problem.b = P_li * problem.b;
    problem.b.conservativeResize(num_li);

    // Remove the components in y corresponding to linearly dependent constraints
    state.y = P_li * state.y;
    state.y.conservativeResize(num_li);

    // Update the names of the dual components y
    if(options.output.active)
		options.output.ynames = extract(options.output.ynames, ili_constraints);
}

auto Regularizer::Impl::fixInfeasibleConstraints(OptimumProblem& problem) -> void
{
    // Auxiliary variables
    const Index m = problem.A.rows();

    // Auxiliary references
    const auto& A = problem.A;
    const auto& l = problem.l;
    auto& b = problem.b;

    // Fix any right-hand side that is infeasible
    for(Index i = 0; i < m; ++i)
        if(min(A.row(i)) >= 0 && min(l) >= 0)
            b[i] = std::max(b[i], dot(A.row(i), l));
        else if(max(A.row(i)) <= 0 && max(l) >= 0)
            b[i] = std::min(b[i], dot(A.row(i), l));
}

auto Regularizer::Impl::echelonizeConstraints(
	OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
{
	// Skip if echelonization should not be performed or A(echelon) was not computed
	if(!params.echelonize || !A_echelon.size())
		return;

	// Compute the corresponding echelonized right-hand side vector b
    problem.b = R * problem.b;

    // Update the y-Lagrange multipliers that correspond now to basic variables
    state.y = tr(invR) * state.y;

    // Update the names of the constraints to the names of basic variables
    if(options.output.active)
        options.output.ynames = extract(options.output.xnames, ibasic_variables);
}

auto Regularizer::Impl::updateConstraints(OptimumProblem& problem) -> void
{
	problem.A = params.echelonize && A_echelon.size() ? A_echelon : A_star;
}

auto Regularizer::Impl::regularize(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
{
	determineTrivialConstraints(problem);
	determineLinearlyDependentConstraints(problem);
	assembleEchelonConstraints(state);

	removeTrivialConstraints(problem, state, options);
	removeLinearlyDependentConstraints(problem, state, options);
	fixInfeasibleConstraints(problem);
	echelonizeConstraints(problem, state, options);
	updateConstraints(problem);
}

auto Regularizer::Impl::recover(const OptimumProblem& problem, OptimumState& state) -> void
{
	// Calculate dual variables y w.r.t. original equality constraints
	if(problem.objective)
		state.y = lu.trsolve(state.f.grad - state.z);
	if(problem.c.size())
		state.y = lu.trsolve(problem.c - state.z);

	// Check if there was any trivial variables and update state accordingly
	if(itrivial_variables.size())
	{
		// Define some auxiliary size variables
		const Index nn = inontrivial_variables.size();
		const Index nt = itrivial_variables.size();
		const Index mn = inontrivial_constraints.size();
		const Index mt = itrivial_constraints.size();
		const Index n = nn + nt;
		const Index m = mn + mt;

		// Set the final objective evaluation
		state.f = f;

		// Resize back the lengths of x, y, z
		state.x.conservativeResize(n);
		state.y.conservativeResize(m);
		state.z.conservativeResize(n);

		// Set the components corresponding to non-trivial variables and constraints
		rows(state.x, inontrivial_variables)   = state.x.segment(0, nn);
		rows(state.y, inontrivial_constraints) = state.y.segment(0, mn);
		rows(state.z, inontrivial_variables)   = state.z.segment(0, nn);

		// Set the components corresponding to trivial variables and constraints
		rows(state.x, itrivial_variables)   = rows(problem.l, itrivial_variables);
		rows(state.y, itrivial_constraints) = 0.0;
		rows(state.z, itrivial_variables)   = 0.0;
	}
}

Regularizer::Regularizer()
: pimpl(new Impl())
{}

Regularizer::Regularizer(const Regularizer& other)
: pimpl(new Impl(*other.pimpl))
{}

Regularizer::~Regularizer()
{}

auto Regularizer::operator=(Regularizer other) -> Regularizer&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Regularizer::regularize(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
{
	pimpl->regularize(problem, state, options);
}

auto Regularizer::recover(const OptimumProblem& problem, OptimumState& state) -> void
{
	pimpl->recover(problem, state);
}

} // namespace Reaktoro
