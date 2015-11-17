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

#include "OptimumSolver.hpp"

// C++ includes
#include <algorithm>
#include <iostream> //todo remove

// Eigen includes
#include <Reaktoro/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumMethod.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolverActNewton.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpAction.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpActive.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpFeasible.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpNewton.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpOpt.hpp>
#include <Reaktoro/Optimization/OptimumSolverKarpov.hpp>
#include <Reaktoro/Optimization/OptimumSolverRefiner.hpp>
#include <Reaktoro/Optimization/OptimumSolverSimplex.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>

namespace Reaktoro {

struct OptimumSolver::Impl
{
    /// The pointer to the optimization solver
    OptimumSolverBase* solver = nullptr;

    /// The IpFeasible solver for approximation calculation
    OptimumSolverIpFeasible ipfeasible;

    // The regularized optimization problem (i.e., no linearly dependent equality constraints and removed trivial constraints)
    OptimumProblem rproblem;

    // The regularized optimization state (i.e., the solution of the regularized optimization problem)
    OptimumState rstate;

    // The regularized optimization options
    OptimumOptions roptions;

    /// The auxiliary objective function evaluations
    ObjectiveResult f, fX;

    // The indices of the linearly independent constraints
    Indices ili_constraints;

    // The indices of the basic variables
    Indices ibasic_variables;

    // The indices of the variables fixed at the lower bound
    Indices itrivial_variables;

    // The indices of the equality constraints whose participating variables are fixed at the lower bound
    Indices itrivial_constraints;

    // The indices of the non-trivial variables
    Indices inontrivial_variables;

    // The indices of the non-trivial constraints
    Indices inontrivial_constraints;

    // The diagonal matrix used to scale the columns of the coefficient matrix `A`
    Vector X;

    // The regularizer matrix that is applied to the original coefficient matrix `A` as `reg(A) = R*A`
    Matrix R;

    // The LU decomposition of the transpose of the coefficient matrix `A`
    Eigen::FullPivLU<Matrix> lu;

    // The lower and upper matrices of the LU decomposition of the coefficient matrix `A`, where `PAQ = LU`
    Matrix L, U;

    // The permutation matrices of the LU decomposition of the coefficient matrix `A`, where `PAQ = LU`
    PermutationMatrix P, Q;

    // Construct a default Impl instance
    Impl()
    {
        // Set the IpNewton method as default
        setMethod(OptimumMethod::IpNewton);
    }

    // Construct a Impl instance with given method
    Impl(OptimumMethod method)
    {
        setMethod(method);
    }

    ~Impl()
    {
        if(solver != nullptr) delete solver;
    }

    // Set the optimization method for the solver
    auto setMethod(OptimumMethod method) -> void
    {
        if(solver != nullptr) delete solver;

        switch(method)
        {
        case OptimumMethod::ActNewton:
            solver = new OptimumSolverActNewton();
            break;
        case OptimumMethod::IpAction:
            solver = new OptimumSolverIpAction();
            break;
        case OptimumMethod::IpActive:
            solver = new OptimumSolverIpActive();
            break;
        case OptimumMethod::IpNewton:
            solver = new OptimumSolverIpNewton();
            break;
        case OptimumMethod::IpOpt:
            solver = new OptimumSolverIpOpt();
            break;
        case OptimumMethod::Karpov:
            solver = new OptimumSolverKarpov();
            break;
        case OptimumMethod::Refiner:
            solver = new OptimumSolverRefiner();
            break;
        case OptimumMethod::Simplex:
            solver = new OptimumSolverSimplex();
            break;
        default:
            solver = new OptimumSolverIpNewton();
            break;
        }
    }

    // Remove trivial constraints and trivial variables  (optional strategy)
    auto regularize1stlevel(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
    {
        // The number of rows and cols in the original coefficient matrix
        const Index m = problem.A.rows();
        const Index n = problem.A.cols();

        // Auxiliary references
        const auto& A = problem.A;
        const auto& b = problem.b;
        const auto& l = problem.l;

        // Return true if the i-th constraint forces the variables to be fixed on the lower bounds
        auto istrivial = [&](Index irow)
        {
            return ( min(A.row(irow)) >= 0 && A.row(irow)*l >= b[irow] ) ||
                   ( max(A.row(irow)) <= 0 && A.row(irow)*l <= b[irow] );
        };

        // Determine the original equality constraints that fix variables on the lower bound
        itrivial_constraints.clear();
        for(Index i = 0; i < m; ++i)
            if(istrivial(i))
                itrivial_constraints.push_back(i);

        // Skip the rest if there is no trivial original constraints
        if(itrivial_constraints.empty())
            return;

        // Determine the trivial original variables that are fixed at their lower bounds
        itrivial_variables.clear();
        for(Index i : itrivial_constraints)
            for(Index j = 0; j < n; ++j)
                if(A(i, j) != 0.0)
                    itrivial_variables.push_back(j);

        // Initialize the indices of the non-trivial original constraints
        inontrivial_constraints = difference(range(m), itrivial_constraints);

        // Initialize the indices of the non-trivial original variables
        inontrivial_variables = difference(range(n), itrivial_variables);

        // Keep only the non-trivial constraints and variables of the original equality constraints
        problem.A = submatrix(problem.A, inontrivial_constraints, inontrivial_variables);
        problem.b = rows(problem.b, inontrivial_constraints);

        // Remove trivial components from problem.objective
        if(problem.objective)
        {
            Vector x = problem.l;

            problem.objective = [=](const Vector& X) mutable
            {
                rows(x, inontrivial_variables) = X;

                f = problem.objective(x);

                fX.val = f.val;
                fX.grad = rows(f.grad, inontrivial_variables);
                fX.hessian.mode = f.hessian.mode;
                if(f.hessian.dense.size())
                    fX.hessian.dense = submatrix(f.hessian.dense, inontrivial_variables, inontrivial_variables);
                if(f.hessian.diagonal.size())
                    fX.hessian.diagonal = rows(f.hessian.diagonal, inontrivial_variables);
                if(f.hessian.inverse.size())
                    fX.hessian.inverse = submatrix(f.hessian.inverse, inontrivial_variables, inontrivial_variables);

                return fX;
            };
        }

        // Remove trivial components from problem.c
        if(problem.c.rows())
            problem.c = rows(problem.c, inontrivial_variables);

        // Remove trivial components from problem.l
        if(problem.l.rows())
            problem.l = rows(problem.l, inontrivial_variables);

        // Remove trivial components from problem.u
        if(problem.u.rows())
            problem.u = rows(problem.u, inontrivial_variables);

        // Keep only non-trivial components corresponding to non-trivial variables and non-trivial constraints
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

    // Remove linearly dependent constraints (non-optional strategy)
    auto regularize2ndlevel(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
    {
        // The transpose of the original coefficient matrix
        Matrix At = tr(problem.A);

        // Initialize the scaling vector `X`
        X = abs(state.x);

        // The threshold used to avoid scaling by very tiny components (prevent it from being zero)
        const double threshold = 1e-10 * (max(X) + 1);

        // Remove very tiny values in X
        X = (X.array() > threshold).select(X, threshold);

        // Scale the columns of matrix A by X
        At = diag(X) * At;

        // Compute the LU decomposition of A with full pivoting (partial piv is not what we want)
        lu = At.fullPivLu();

        // The rank of the 1st level regularized coefficient matrix A
        const Index rank = lu.rank();

        // Initialize the L, U, P, Q matrices
        L = tr(lu.matrixLU()).topLeftCorner(rank, rank).triangularView<Eigen::Lower>();
        U = tr(lu.matrixLU()).topRows(rank).triangularView<Eigen::UnitUpper>();
        P = lu.permutationQ().inverse();
        Q = lu.permutationP().inverse();

        // Correct the U matrix by unscaling it by X
        U = U * Q.inverse() * diag(inv(X)) * Q;

        std::cout << "PAQ - LU = \n" << P*problem.A*Q - L*U << std::endl;

        // Initialize the indices of the original constraints that are linearly independent
        ili_constraints = Indices(P.indices().data(), P.indices().data() + rank);
        
//        // Sort the indices of the linearly independent constraints
//        std::sort(ili_constraints.begin(), ili_constraints.end());

        // Initialize the indices of the basic variables
        ibasic_variables = Indices(Q.indices().data(), Q.indices().data() + rank);

        // Remove linearly dependent rows from A and b
        problem.A = P * problem.A;
        problem.A = rows(problem.A, 0, rank);

        problem.b = P * problem.b;
        problem.b = rows(problem.b, 0, rank);

//        problem.A = rows(problem.A, ili_constraints);
//        problem.b = rows(problem.b, ili_constraints);

        // Keep only components that correspond to linearly independent constraints
        state.y = rows(state.y, ili_constraints);
    }

    // Transform the equality constraints into cannonical form (optional strategy)
    auto regularize3rdlevel(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void
    {
        // Auxiliary variables
        const Index m = problem.A.rows();
        const Indices& B = ibasic_variables;

        // Update the y-Lagrange multipliers to the residuals of the basic variables (before A is changed!)
        state.y.resize(m);
        for(Index i = 0; i < m; ++i)
            state.y[i] = tr(problem.A).row(B[i]) * state.y;

        // The rank of the 1st level regularized coefficient matrix A
        const Index rank = lu.rank();

        // Create a reference to the U1 part of U = [U1 U2]
        const auto U1 = U.leftCols(rank).triangularView<Eigen::Upper>();

        // Compute the regularizer matrix R
        R = L.triangularView<Eigen::Lower>().solve(identity(rank, rank));
        R = U1.solve(R);

        // Compute the 2nd level equality constraint regularization
        problem.A = R * problem.A;
        problem.b = R * problem.b;

        // Check if the regularizer matrix is composed of
        // rationals that can be recovered from round-off errors
        if(options.max_denominator)
        {
            cleanRationalNumbers(problem.A, options.max_denominator);
            cleanRationalNumbers(R, options.max_denominator);
        }

        // Update the names of the constraints
        if(options.output.active)
            options.output.ynames = extract(options.output.xnames, ibasic_variables);
    }

    // Ensure no positive or negative constraints have infeasible right-hand side
    auto regularize4thlevel(OptimumProblem& problem) -> void
    {
        // Auxiliary variables
        const Index m = problem.A.rows();

        // Auxiliary references
        const auto& A = problem.A;
        const auto& b = problem.b;
        const auto& l = problem.l;

        std::cout << "A = \n" << A << std::endl;
        std::cout << "b = \n" << b << std::endl;

        // Fix any right-hand side that is infeasible
        for(Index i = 0; i < m; ++i)
            if(min(A.row(i)) >= 0 && min(l) >= 0)
                problem.b[i] = std::max(b[i], dot(A.row(i), l));
            else if(max(A.row(i)) <= 0 && max(l) >= 0)
                problem.b[i] = std::min(b[i], dot(A.row(i), l));
    }

    auto initialize(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> void
    {
        // Assert either objective function or c vector have been initialized
        if(problem.objective && problem.c.size())
            RuntimeError("Cannot solve the optimization problem.",
                "The given OptimumProblem instance is ambiguous: "
                    "both members `c` and `objective` have been initialized.");

        // Auxiliary variables
        const auto m = problem.A.rows();
        const auto n = problem.A.cols();

        // Ensure x, y, z have correct dimensions
        if(state.x.size() != n) state.x = zeros(n);
        if(state.y.size() != m) state.y = zeros(m);
        if(state.z.size() != n) state.z = zeros(n);

        // Initialize the regularized problem, state, and options instances
        rproblem = problem;
        rstate = state;
        roptions = options;

        regularize1stlevel(rproblem, rstate, roptions);
        regularize2ndlevel(rproblem, rstate, roptions);
        regularize3rdlevel(rproblem, rstate, roptions);
        regularize4thlevel(rproblem);
    }

    auto finalize(const OptimumProblem& problem, OptimumState& state) -> void
    {
        if(problem.objective)
            rstate.y = lu.solve(rstate.f.grad - rstate.z);

        if(problem.c.size())
            rstate.y = lu.solve(rproblem.c - rstate.z);

        if(itrivial_variables.empty())
        {
            state = rstate;
        }
        else
        {
            state.f = f;

            rows(state.x, inontrivial_variables) = rstate.x;
            rows(state.y, inontrivial_constraints) = rstate.y;
            rows(state.z, inontrivial_variables) = rstate.z;

            rows(state.x, itrivial_variables) = rows(problem.l, itrivial_variables);
            rows(state.y, itrivial_constraints) = 0.0;
            rows(state.z, itrivial_variables) = 0.0;
        }

    }

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        initialize(problem, state, options);

        OptimumResult res = solver->solve(rproblem, rstate, roptions);

        finalize(problem, state);

        return res;
    }
};

OptimumSolver::OptimumSolver()
: pimpl(new Impl())
{}

OptimumSolver::OptimumSolver(OptimumMethod method)
: pimpl(new Impl(method))
{}

OptimumSolver::OptimumSolver(const OptimumSolver& other)
: pimpl(new Impl(*other.pimpl))
{
    // Clone the optimization solver from the copying object
    pimpl->solver = other.pimpl->solver->clone();
}

OptimumSolver::~OptimumSolver()
{}

auto OptimumSolver::operator=(OptimumSolver other) -> OptimumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolver::setMethod(OptimumMethod method) -> void
{
    pimpl->setMethod(method);
}

auto OptimumSolver::approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return approximate(problem, state, {});
}

auto OptimumSolver::approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->ipfeasible.approximate(problem, state, options);
}

auto OptimumSolver::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return solve(problem, state, {});
}

auto OptimumSolver::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
