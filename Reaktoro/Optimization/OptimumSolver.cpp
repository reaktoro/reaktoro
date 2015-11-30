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

// Eigen includes
#include <Reaktoro/Eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Math/LU.hpp>
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

    // The auxiliary objective function evaluations
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

    // The original coefficient matrix `A` used last time
    Matrix A;

    // The diagonal matrix used to scale the columns of the original coefficient matrix `A`
    Vector X;

    // The regularizer matrix that is applied to the coefficient matrix `Abar` with linearly
    // independent rows and non-trivial constraints and variables as `Areg = R*Abar`
    Matrix R;

    // The inverse of the regularizer matrix R
    Matrix invR;

    // The full-pivoting LU decomposition of the coefficient matrix A
    LU lu;

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
    auto regularize1stlevel(const OptimumProblem& problem) -> void
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

        // Skip the rest if there is no trivial constraints
        if(itrivial_constraints.empty())
            return;

        // Determine the trivial original variables that are fixed at their lower bounds
        itrivial_variables.clear();
        for(Index i = 0; i < n; ++i)
            for(Index j = 0; j < m; ++j)
                if(A(j, i) != 0.0 && contained(j, itrivial_constraints))
                    { itrivial_variables.push_back(i); break; }

        // Initialize the indices of the non-trivial original constraints
        inontrivial_constraints = difference(range(m), itrivial_constraints);

        // Initialize the indices of the non-trivial original variables
        inontrivial_variables = difference(range(n), itrivial_variables);

        // Assert there not all contraints are trivial
        Assert(inontrivial_variables.size(),
            "Could not perform the optimization calculation.",
            "The provided optimization problem contains only trivial constraints.");

        // Remove trivial constraints and variables of the equality constraints
        rproblem.A = submatrix(problem.A, inontrivial_constraints, inontrivial_variables);
        rproblem.b = rows(problem.b, inontrivial_constraints);

        // Remove trivial components from problem.c
        if(problem.c.rows())
            rproblem.c = rows(problem.c, inontrivial_variables);

        // Remove trivial components from problem.l
        if(problem.l.rows())
            rproblem.l = rows(problem.l, inontrivial_variables);

        // Remove trivial components from problem.u
        if(problem.u.rows())
            rproblem.u = rows(problem.u, inontrivial_variables);

        // The auxiliary vector used in the lambda functions below. The use of this vector `x`
        // ensures that trivial components remain unchanged on the lower bounds.
        Vector x = problem.l;

        // Remove trivial components from problem.objective
        if(problem.objective)
            rproblem.objective = [=](const Vector& X) mutable
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

        // Remove trivial components corresponding to trivial variables and trivial constraints
        rstate.x = rows(rstate.x, inontrivial_variables);
        rstate.y = rows(rstate.y, inontrivial_constraints);
        rstate.z = rows(rstate.z, inontrivial_variables);

        // Update the names of the constraints and variables accordingly
        if(roptions.output.active)
        {
            roptions.output.xnames = extract(roptions.output.xnames, inontrivial_variables);
            roptions.output.ynames = extract(roptions.output.ynames, inontrivial_constraints);
            roptions.output.znames = extract(roptions.output.znames, inontrivial_variables);
        }
    }

    // Remove linearly dependent constraints (non-optional strategy)
    auto regularize2ndlevel() -> void
    {
        // Auxiliary variables
        const Index n = rproblem.A.cols();

        // Initialize the scaling vector `X`
        X = abs(rstate.x);

        // The threshold used to avoid scaling by very tiny components (prevent it from being zero)
        const double threshold = 1e-10 * (max(X) + 1);

        // Remove very tiny values in X
        X = (X.array() > threshold).select(X, threshold);

        // Compute the LU decomposition of coefficient matrix A
        lu.compute(rproblem.A, X);

        // Auxiliary references to LU components
        const auto& P = lu.P;
        const auto& Q = lu.Q;
        const auto& rank = lu.rank;

        // Initialize the indices of the original constraints that are linearly independent
        ili_constraints = Indices(P.indices().data(), P.indices().data() + rank);
        
        // Initialize the indices of the basic variables
        ibasic_variables = Indices(Q.indices().data(), Q.indices().data() + rank);

        // Permute the rows of A and b
        rproblem.A = P * rproblem.A;
        rproblem.b = P * rproblem.b;

        // Remove the rows of A and b past rank
        rproblem.A.conservativeResize(rank, n);
        rproblem.b.conservativeResize(rank);

        // Keep only components that correspond to linearly independent constraints
        rstate.y = P * rstate.y;
        rstate.y.conservativeResize(rank);
    }

    // Transform the equality constraints into cannonical form (optional strategy)
    auto regularize3rdlevel(const OptimumProblem& problem, const OptimumState& state) -> void
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

        // Compute the 2nd level equality constraint regularization
        rproblem.A = R * rproblem.A;

        // Check if the regularizer matrix is composed of rationals that can be recovered from round-off errors
        if(roptions.max_denominator)
        {
            cleanRationalNumbers(rproblem.A, roptions.max_denominator);
            cleanRationalNumbers(R, roptions.max_denominator);
            cleanRationalNumbers(invR, roptions.max_denominator);
        }

        // Update the y-Lagrange multipliers to the residuals of the basic variables (before A is changed!)
        rstate.y = tr(invR) * rstate.y;

        // After round-off cleanup, compute the 2nd level regularization of b
        rproblem.b = R * rproblem.b;

        // Update the names of the constraints
        if(roptions.output.active)
            roptions.output.ynames = extract(roptions.output.xnames, ibasic_variables);
    }

    // Ensure no positive or negative constraints have infeasible right-hand side
    auto regularize4thlevel() -> void
    {
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

    // Initialize the regularized rproblem, rstate, and roptions objects
    auto initialize(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> void
    {
        // Assert either objective function or c vector have been initialized
        if(problem.objective && problem.c.size())
            RuntimeError("Cannot solve the optimization problem.",
                "The given OptimumProblem instance is ambiguous: "
                    "both members `c` and `objective` have been initialized.");

        // The number of variables, constraints, and sensitivity parameters
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

        // Regularize the equality constraints in four levels
        regularize1stlevel(problem);
        regularize2ndlevel();
        regularize3rdlevel(problem, state);
        regularize4thlevel();
    }

    // Trasnfer the solution state of the regularized rstate to state
    auto finalize(const OptimumProblem& problem, OptimumState& state) -> void
    {
        // Calculate dual variables y w.r.t. original equality constraints
        if(problem.objective)
            rstate.y = lu.trsolve(rstate.f.grad - rstate.z);
        if(problem.c.size())
            rstate.y = lu.trsolve(rproblem.c - rstate.z);

        // Check if there was any trivial variables and update state accordingly
        if(itrivial_variables.empty())
        {
            state = rstate;
        }
        else
        {
            // Set the final objective evaluation
            state.f = f;

            // Set the components corresponding to non-trivial variables and constraints
            rows(state.x, inontrivial_variables) = rstate.x;
            rows(state.y, inontrivial_constraints) = rstate.y;
            rows(state.z, inontrivial_variables) = rstate.z;

            // Set the components corresponding to trivial variables and constraints
            rows(state.x, itrivial_variables) = rows(problem.l, itrivial_variables);
            rows(state.y, itrivial_constraints) = 0.0;
            rows(state.z, itrivial_variables) = 0.0;
        }
    }

    // Solve the optimization problem
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        initialize(problem, state, options);

        OptimumResult res = solver->solve(rproblem, rstate, roptions);

        finalize(problem, state);

        return res;
    }

    /// Calculate the sensitivity of the optimal solution with respect to parameters.
    auto dxdp(Vector dgdp, Vector dbdp) -> Vector
    {
        // Assert the size of the input matrices dgdp and dbdp
        Assert(dgdp.size() && dbdp.size() && dgdp.cols() == dbdp.cols(),
            "Could not calculate the sensitivity of the optimal solution with respect to parameters.",
            "The given input matrices `dgdp` and `dbdp` are either empty or does not have the same number of columns.");

        const Index n = dgdp.rows();

        const auto& r = lu.rank;
        const auto& P = lu.P;

        // Remove derivative components corresponding to trivial constraints
        if(itrivial_constraints.size())
        {
            dbdp = rows(dbdp, inontrivial_constraints);
            dgdp = rows(dgdp, inontrivial_variables);
        }

        // Permute the rows of dbdp and remove its rows past rank
        dbdp = P * dbdp;
        dbdp.conservativeResize(r);

        // After round-off cleanup, compute the 2nd level regularization of dbdp
        dbdp = R * dbdp;

        // Set the components of the sensitivity matrix corresponding to trivial and non-trivial variables
        if(itrivial_constraints.size())
        {
            Vector dxdp(n);
            auto aux = solver->dxdp(dgdp, dbdp);
            rows(dxdp, inontrivial_variables) = solver->dxdp(dgdp, dbdp);
            rows(dxdp, itrivial_variables) = 0.0;
            return dxdp;
        }
        else return solver->dxdp(dgdp, dbdp);
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

auto OptimumSolver::dxdp(const Vector& dgdp, const Vector& dbdp) -> Vector
{
    return pimpl->dxdp(dgdp, dbdp);
}

} // namespace Reaktoro
