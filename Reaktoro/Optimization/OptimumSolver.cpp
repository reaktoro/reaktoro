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

    // The QR decomposition of the transpose of coefficient matrix `A`
    Eigen::ColPivHouseholderQR<Matrix> qr;

    // The regularized optimization problem (i.e., no linearly dependent equality constraints and removed trivial constraints)
    OptimumProblem regularized_problem;

    // The regularized optimization state (i.e., the solution of the regularized optimization problem)
    OptimumState regularized_state;

    // The regularized optimization options
    OptimumOptions regularized_options;

    /// The auxiliary objective function evaluations
    ObjectiveResult f, fX;

    // The indices of the working variables
    Indices iworking_variables;

    // The indices of the working constraints
    Indices iworking_constraints;

    /// The indices of the equality constraints with positive only coefficients
    Indices iconstraints_positive_coefficients;

    /// The indices of the equality constraints with negative only coefficients
    Indices iconstraints_negative_coefficients;

    // The indices of the variables fixed at the lower bound
    Indices itrivial_variables;

    // The indices of the equality constraints whose participating variables are fixed at the lower bound
    Indices itrivial_constraints;

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

    // Determine the equality constraints that always require positive right-hand sides
    auto determinePositiveAndNegativeConstraints(const OptimumProblem& problem) -> void
    {
        const Matrix& A = problem.A;
        iconstraints_positive_coefficients.clear();
        iconstraints_negative_coefficients.clear();
        for(int i = 0; i < A.rows(); ++i)
            if(min(A.row(i)) >= 0)
                iconstraints_positive_coefficients.push_back(i);
            else if(max(A.row(i)) <= 0)
                iconstraints_negative_coefficients.push_back(i);
            else continue;
    }

    // Determine the trivial variables that are fixed at their lower bounds
    auto determineTrivialVariablesAndConstraints(const OptimumProblem& problem) -> void
    {
        const Matrix& A = problem.A;
        const Vector& b = problem.b;
        const Vector& l = problem.l;
        itrivial_variables.clear();
        itrivial_constraints.clear();
        for(Index i : iconstraints_positive_coefficients)
            if(A.row(i)*l >= b[i])
                itrivial_constraints.push_back(i);
        for(Index i : iconstraints_negative_coefficients)
            if(A.row(i)*l <= b[i])
                itrivial_constraints.push_back(i);
        for(Index i : itrivial_constraints)
            for(int j = 0; j < A.cols(); ++j)
                if(A(i, j) != 0.0)
                    itrivial_variables.push_back(j);
    }

    auto determineWorkingVariablesAndConstraints(const OptimumProblem& problem) -> void
    {
        const Matrix& A = problem.A;

        // Compute the QR factorization of the transpose of the coefficient matrix `A`
        qr.compute(tr(A));

        // Identify the indices of the linearly independent rows of `A`
        const unsigned rank = qr.rank();
        Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
        std::sort(I.data(), I.data() + rank);

        // The indices of the linearly independent rows of `A`
        iworking_constraints = Indices(I.data(), I.data() + rank);

        // Remove the trivial equality constraints (i.e., those where the participating variables are hold at their lower bounds)
        iworking_constraints = difference(iworking_constraints, itrivial_constraints);

        // The indices of the working variables (i.e., all variables except the trivial ones)
        iworking_variables = difference(range<Index>(A.cols()), itrivial_variables);
    }

    auto regularizeOptimumProblem(const OptimumProblem& problem) -> void
    {
        if(problem.objective && problem.c.size())
            RuntimeError("Cannot solve the optimization problem.",
                "The given OptimumProblem instance is ambiguous: "
                    "both members `c` and `objective` have been defined.");
        determinePositiveAndNegativeConstraints(problem);
        determineTrivialVariablesAndConstraints(problem);
        determineWorkingVariablesAndConstraints(problem);
        regularized_problem.A = submatrix(problem.A, iworking_constraints, iworking_variables);
        regularized_problem.b = rows(problem.b, iworking_constraints);
        if(problem.c.rows())
            regularized_problem.c = rows(problem.c, iworking_variables);
        if(problem.l.rows())
            regularized_problem.l = rows(problem.l, iworking_variables);
        if(problem.u.rows())
            regularized_problem.u = rows(problem.u, iworking_variables);
        if(problem.objective)
        {
            Vector x = problem.l;

            regularized_problem.objective = [=](const Vector& X) mutable
            {
                rows(x, iworking_variables) = X;

                f = problem.objective(x);

                fX.val = f.val;
                fX.grad = rows(f.grad, iworking_variables);
                fX.hessian.mode = f.hessian.mode;
                if(f.hessian.dense.size())
                    fX.hessian.dense = submatrix(f.hessian.dense, iworking_variables, iworking_variables);
                if(f.hessian.diagonal.size())
                    fX.hessian.diagonal = rows(f.hessian.diagonal, iworking_variables);
                if(f.hessian.inverse.size())
                    fX.hessian.inverse = submatrix(f.hessian.inverse, iworking_variables, iworking_variables);

                return fX;
            };
        }
    }

    auto regularizeOptimumOptions(const OptimumOptions& options) -> void
    {
        regularized_options = options;

        // Keep only the names of the working variables
        if(options.output.xnames.size())
            regularized_options.output.xnames = extract(options.output.xnames, iworking_variables);

        // Keep only the names of the working variables
        if(options.output.znames.size())
            regularized_options.output.znames = extract(options.output.znames, iworking_variables);

        // Keep only the names of the working constraints
        if(options.output.ynames.size())
            regularized_options.output.ynames = extract(options.output.ynames, iworking_constraints);
    }

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        regularizeOptimumProblem(problem);
        regularizeOptimumOptions(options);

        state.x.resize(problem.A.cols());
        state.y.resize(problem.A.rows());
        state.z.resize(problem.A.cols());

        regularized_state.x = rows(state.x, iworking_variables);
        regularized_state.y = rows(state.y, iworking_constraints);
        regularized_state.z = rows(state.z, iworking_variables);

        OptimumResult res = solver->solve(regularized_problem, regularized_state, regularized_options);

        rows(state.x, iworking_variables) = regularized_state.x;
        rows(state.y, iworking_constraints) = regularized_state.y;
        rows(state.z, iworking_variables) = regularized_state.z;

        state.f = f;

        rows(state.x, itrivial_variables) = rows(problem.l, itrivial_variables);
        rows(state.y, itrivial_constraints) = 0.0;
        rows(state.z, itrivial_variables) = 0.0;

        if(problem.objective)
            state.y = qr.solve(f.grad - state.z);

        if(problem.c.size())
            state.y = qr.solve(problem.c - state.z);

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
