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

#include "RegularizationUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>

namespace Reaktoro {

auto Regularizer::regularizeTrivialConstraints() -> void
{
    // Auxiliary references
    const Matrix& A = problem.A;
    const Vector& b = problem.b;
    const Vector& l = problem.l;

    // The number of rows and cols in the original coefficient matrix
    const Index m = A.rows();
    const Index n = A.cols();

    // Only determine trivial constraints if A is different than last A
    if(problem.A != A_last)
    {
        // Return true if the i-th constraint forces the variables to be fixed on the lower bounds
        auto istrivial = [&](Index irow)
        {
            return ( min(A.row(irow)) >= 0 && A.row(irow)*l >= b[irow] ) ||
                   ( max(A.row(irow)) <= 0 && A.row(irow)*l <= b[irow] );
        };

        // Clear previous states of trivial and non-trivial constraints and variables
        itrivial_constraints.clear();
        itrivial_variables.clear();
        inontrivial_constraints.clear();
        inontrivial_variables.clear();

        // Determine the original equality constraints that fix variables on the lower bound
        for(Index i = 0; i < m; ++i)
            if(istrivial(i))
                itrivial_constraints.push_back(i);

        // Check if there are no trivial constraints
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

            // Update the matrix `A_star`
            A_star = submatrix(A, inontrivial_constraints, inontrivial_variables);

            // Assert there not all contraints are trivial
            Assert(inontrivial_variables.size(),
                "Could not accept the optimization problem.",
                "The provided problem contains only trivial constraints.");
        }
        else
        {
            A_star = A;
        }
    }

    // Check if there is any trivial constraint. If yes, then update problem members.
    if(itrivial_constraints.empty())
    {
        // no need to update A_star, since it has not changed or was updated above
        b_star = b;
    }
    else
    {
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

        // Update the b_star (no need to update A_star as it did not change or was updated above)
        b_star = rows(b, inontrivial_constraints);

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
}

auto Regularizer::regularizeLinearlyDependentConstraints() -> void
{
    // Auxiliary references to LU components
    const auto& P = lu.P;
    const auto& Q = lu.Q;
    const auto& rank = lu.rank;

    // Only determine linearly dependent rows if A is different than last A
    if(problem.A == A_last)
    {
        // Keep only components that correspond to linearly independent constraints
        state.y = P * state.y;
        state.y.conservativeResize(rank);
        return;
    }

    // The number of rows and cols in the original coefficient matrix
    const Index m = A_star.rows();
    const Index n = A_star.cols();

    // Compute the LU decomposition of coefficient matrix A
    lu.compute(A_star);

    // Check if all rows are linearly independent
    if(rank != m)
    {
        // Update the permutation matrix and the number of linearly independent constraints
        P_li = lu.P;
        num_li = lu.rank;

        // Permute the rows of A and b
        A_star = P_li * A_star;
        b_star = P_li * b_star;

        // Remove the rows of A and b past rank
        A_star.conservativeResize(num_li, n);
        b_star.conservativeResize(num_li);

        // Keep only components that correspond to linearly independent constraints
        state.y = P_li * state.y;
        state.y.conservativeResize(num_li);
    }
    else
    {
        // no need to update b_star and A_star here, since it has all rows in A_star are linearly independent
        P_li.setIdentity(m);
        num_li = m;
    }
}

auto Regularizer::regularizeNonCannonicalConstraints() -> void
{
    // Auxiliary variables
    const Index n = A_star.cols();

    // Initialize the weight vector `W`
    W = abs(state.x);

    // The threshold used to avoid scaling by very tiny components (prevent it from being zero)
    const double threshold = 1e-10 * (max(W) + 1);

    // Remove very tiny values in W
    W = (W.array() > threshold).select(W, threshold);

    // Compute the LU decomposition of matrix A_star with column-sorting weights.
    // Columsn corresponding to variables with higher weights are moved to the beginning of the matrix.
    // This gives preference for those variables to become linearly independent basis
    lu.compute(A_star, W);

    // Auxiliary references to LU components
    const auto& P = lu.P;
    const auto& Q = lu.Q;
    const auto& rank = lu.rank;

    // Initialize the indices of the basic variables
    ibasic_variables = Indices(Q.indices().data(), Q.indices().data() + rank);

    // Check if the new set of basic variables is diffent than the previous
    if(!contained(ibasic_variables, ibasic_variables_last))
    {
        // Update the set of last basic variables
        ibasic_variables_last = ibasic_variables;

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
        A_reg = R * A_star;

        // Check if the regularizer matrix is composed of rationals that can be recovered from round-off errors
        if(options.max_denominator)
        {
            cleanRationalNumbers(A_reg, options.max_denominator);
            cleanRationalNumbers(R, options.max_denominator);
            cleanRationalNumbers(invR, options.max_denominator);
        }

        // Update the names of the constraints
        if(options.output.active)
            options.output.ynames = extract(options.output.xnames, ibasic_variables);

        // Update the A member of problem
        problem.A = A_reg;
    }

    // Update the y-Lagrange multipliers to the residuals of the basic variables (before A is changed!)
    state.y = tr(invR) * state.y;

    // After round-off cleanup, compute the regularization of b
    problem.b = b_reg = R * b_star;
}

auto Regularizer::regularizeInfeasibleConstraints() -> void
{
    auto
    // Auxiliary variables
    const Index m = rproblem.A.rows();

    // Auxiliary references
    const auto& A = rproblem.A;
    const auto& b = rproblem.b;
    const auto& l = rproblem.l;

    // Fix any right-hand side that is infeasible
    for(Index i = 0; i < m; ++i)
        if(min(A.row(i)) >= 0 && min(l) >= 0)
            rproblem.b[i] = std::max(b[i], dot(A.row(i), l));
        else if(max(A.row(i)) <= 0 && max(l) >= 0)
            rproblem.b[i] = std::min(b[i], dot(A.row(i), l));
}

auto RegOptimumProblem::regularizeOptimumState(const OptimumState& orig, OptimumState& reg) -> void
{
    // Remove trivial components corresponding to trivial variables and trivial constraints
    if(inontrivial_constraints.size())
    {
        reg.x = rows(orig.x, inontrivial_variables);
        reg.y = rows(orig.y, inontrivial_constraints);
        reg.z = rows(orig.z, inontrivial_variables);
    }

    // Keep only components that correspond to linearly independent constraints
    reg.y = P * orig.y;
    reg.y.conservativeResize(rank);
}

auto RegOptimumProblem::regularizeOptimumOptions(const OptimumOptions& orig, OptimumOptions& reg) -> void
{

    // Update the names of the constraints and variables accordingly
    if(inontrivial_constraints.size() && orig.output.active)
    {
        reg.output.xnames = extract(orig.output.xnames, inontrivial_variables);
        reg.output.ynames = extract(orig.output.ynames, inontrivial_constraints);
        reg.output.znames = extract(orig.output.znames, inontrivial_variables);
    }
}

} // namespace Reaktoro
