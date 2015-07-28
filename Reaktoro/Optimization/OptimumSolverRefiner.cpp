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

#include "OptimumSolverRefiner.hpp"

// C++ includes
#include <iostream>

// Eigen includes
#include <eigen/Dense>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {

struct OptimumSolverRefiner::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;
    ObjectiveResult f;

    Outputter outputter;

    auto solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult;

    auto solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverRefiner::Impl::solve(OptimumProblem problem, OptimumState& state, OptimumOptions options) -> OptimumResult
{
    // The transpose of the coefficient matrix `A`
    const Matrix At = tr(problem.A);

    // Calculate the QR decomposition of the transpose of `A`
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(At);

    // Identify the indices of the linearly independent rows of `A`
    const unsigned rank = qr.rank();
    Eigen::VectorXi I = qr.colsPermutation().indices().segment(0, rank);
    std::sort(I.data(), I.data() + rank);

    // The indices of the linearly independent rows of `A`
    const Indices ic(I.data(), I.data() + rank);

    // Define the regularized optimization problem without linearly dependent constraints
    problem.A = rows(problem.A, ic);
    problem.b = rows(problem.b, ic);

    // Remove the names of the linearly dependent constraints
    if(options.output.ynames.size())
        options.output.ynames = extract(options.output.ynames, ic);

    // Get the linearly independent components of the Lagrange multipliers `y`
    state.y = rows(state.y, ic);

    // Solve the regularized optimization problem
    auto result = solveMain(problem, state, options);

    // Calculate the Lagrange multipliers for all equality constraints
    state.y = qr.solve(f.grad - state.z);

    return result;
}

auto OptimumSolverRefiner::Impl::solveMain(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Start timing the calculation
    Time begin = time();

    // Initialize the outputter instance
    outputter = Outputter();
    outputter.setOptions(options.output);

    // Set the KKT options
    kkt.setOptions(options.kkt);

    // The result of the calculation
    OptimumResult result;

    // Define some auxiliary references to variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;

    // Define some auxiliary references to parameters
    const auto& n         = problem.A.cols();
    const auto& m         = problem.A.rows();
    const auto& tolerance = options.tolerance;
    const auto& tau       = options.ipnewton.tau;

    auto A = problem.A;
    auto b = problem.b;

    // The transpose representation of matrix `A`
    const auto At = tr(A);

    auto lu = A.fullPivLu();

    // Get the lower and upper matrices
    Matrix L = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitLower>();
    Matrix U = lu.matrixLU().triangularView<Eigen::Upper>();

    // Get the permutation matrices
    const auto P1 = lu.permutationP();
    const auto P2 = lu.permutationQ();

    // Set the U1 and U2 submatrices of U = [U1 U2]
    const Matrix U1 = U.leftCols(m);
    const Matrix U2 = U.rightCols(n - m);

    A = U*P2;
    b = P1.inverse() * b;
    b = L.triangularView<Eigen::UnitLower>().solve(b);

//    Matrix L  = lu.matrixLU().leftCols(m).triangularView<Eigen::Lower>();
//    Matrix U1 = lu.matrixLU().leftCols(m).triangularView<Eigen::UnitUpper>();
//    Matrix U2 = lu.matrixLU().leftCols(n - m);
//
//    Matrix U(m, n); U << U1, U2;
//
//    std::cout << L*U << std::endl;



    Vector h;

    // Ensure the initial guesses for `x` and `y` have adequate dimensions
    if(x.size() != n) x = zeros(n);
    if(y.size() != m) y = zeros(m);

    // Ensure the dual potentials `z` are zero
    z = zeros(n);

    // Ensure the initial guess for `x` is inside the feasible domain
    x = (x.array() > 0.0).select(x, 1.0);

    // Calculate the kernel (nullspace) matrix K
    const Matrix K = A.fullPivLu().kernel();

    // The alpha step size used to restric the steps inside the feasible domain
    double alpha;

    // The optimality, feasibility, and total error variables
    double errorf, errorh, error;

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(!options.output.active) return;

        outputter.addEntry("iter");
        outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        outputter.addEntries(options.output.yprefix, m, options.output.ynames);
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("errorf");
        outputter.addEntry("errorh");
        outputter.addEntry("errorc");
        outputter.addEntry("error");
        outputter.addEntry("alpha");

        outputter.outputHeader();
        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.outputState();
    };

    // The function that outputs the current state of the solution
    auto output_state = [&]()
    {
        if(!options.output.active) return;

        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(error);
        outputter.addValue(alpha);
        outputter.outputState();
    };

    // The function that updates the objective and constraint state
    auto update_state = [&]()
    {
        f = problem.objective(x);
        h = A*x - b;
    };

    // Return true if function `update_state` failed
    auto update_state_failed = [&]()
    {
        const bool f_finite = std::isfinite(f.val);
        const bool g_finite = f.grad.allFinite();
        const bool all_finite = f_finite && g_finite;
        return !all_finite;
    };

    // The function that computes the Newton step
    auto compute_newton_step_gem = [&]()
    {
//        Matrix J = zeros(n, n);
//        block(J, 0, 0, n - m, n) = tr(K) * diag(f.hessian.diagonal);
//        block(J, n - m, 0, m, n) = A;
//
//        Matrix r = zeros(n);
//        rows(r, 0, n - m) = -tr(K) * f.grad;
//        rows(r, n - m, m) = -(A*x - b);
//
//        sol.dx = J.lu().solve(r);
//        sol.dy = tr(A).fullPivLu().solve(f.hessian.diagonal % sol.dx + f.grad - At*y);

        Matrix J = zeros(n + m, n + m);
        block(J, 0, 0, n, n) = diag(f.hessian.diagonal);
        block(J, 0, n, n, m) = -At;
        block(J, n, 0, m, n) =  A;

        Matrix r = zeros(n + m);
        rows(r, 0, n)= -(f.grad - At*y);
        rows(r, n, m)= -h;

        Vector u = J.lu().solve(r);

        sol.dx = rows(u, 0, n);
        sol.dy = rows(u, n, m);

//        KktMatrix lhs{f.hessian, A, x, z};
//
//        kkt.decompose(lhs);
//
//        // Compute the right-hand side vectors of the KKT equation
//        rhs.rx.noalias() = -(f.grad - At*y);
//        rhs.ry.noalias() = -(h);
//        rhs.rz.noalias() = zeros(n);
//
//        // Compute `dx` and `dy` by solving the KKT equation
//        kkt.solve(rhs, sol);
//
//        // Update the time spent in linear systems
//        result.time_linear_systems += kkt.result().time_solve;
//        result.time_linear_systems += kkt.result().time_decompose;
    };

    // The function that computes the Newton step
    auto compute_newton_step_lma = [&]()
    {
        Matrix J = zeros(n, n);
        block(J, 0, 0, n - m, n) = tr(K) * diag(f.hessian.diagonal);
        block(J, n - m, 0, m, n) = A;

        Matrix r = zeros(n);
        rows(r, 0, n - m) = -tr(K) * f.grad;
        rows(r, n - m, m) = -(A*x - b);

        sol.dx = J.lu().solve(r);
    };

    // Return true if the function `compute_newton_step` failed
    auto compute_newton_step_failed_lma = [&]()
    {
        return !sol.dx.allFinite();
    };

    // Return true if the function `compute_newton_step` failed
    auto compute_newton_step_failed_gem = [&]()
    {
        const bool dx_finite = sol.dx.allFinite();
        const bool dy_finite = sol.dy.allFinite();
        const bool all_finite = dx_finite && dy_finite;
        return !all_finite;
    };

    // The function that performs an update in the iterates
    auto update_iterates_lma = [&]()
    {
        alpha = fractionToTheBoundary(x, sol.dx, tau);
        x += alpha * sol.dx;
    };

    // The function that performs an update in the iterates
    auto update_iterates_gem = [&]()
    {
        alpha = fractionToTheBoundary(x, sol.dx, tau);

        if(options.ipnewton.uniform_newton_step)
        {
            x += alpha * sol.dx;
            y += alpha * sol.dy;
        }
        else
        {
            x += alpha * sol.dx;
            y += sol.dy;
        }
    };

    // The function that computes the current error norms
    auto update_errors_gem = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
        errorf = norminf(f.grad - At*y);
        errorh = norminf(h);

        // Calculate the maximum error
        error = std::max(errorf, errorh);
        result.error = error;
    };

    // The function that computes the current error norms
    auto update_errors_lma = [&]()
    {
        // Calculate the optimality and feasibility errors
        errorf = norminf(tr(K) * f.grad);
        errorh = norminf(h);

        // Calculate the maximum error
        error = std::max(errorf, errorh);
        result.error = error;
    };

    auto converged = [&]()
    {
        if(error < tolerance)
        {
            result.succeeded = true;
            return true;
        }
        return false;
    };

    update_state();
    output_header();

    do
    {
        ++result.iterations; if(result.iterations > options.max_iterations) break;

        if(options.refiner.use_lma_setup)
        {
            compute_newton_step_lma();
            if(compute_newton_step_failed_lma())
                break;
            update_iterates_lma();
            update_state();
            if(update_state_failed())
                break;
            update_errors_lma();
        }
        else
        {
            compute_newton_step_gem();
            if(compute_newton_step_failed_gem())
                break;
            update_iterates_gem();
            update_state();
            if(update_state_failed())
                break;
            update_errors_gem();
        }
            output_state();
    } while(!converged());

    outputter.outputHeader();

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverRefiner::OptimumSolverRefiner()
: pimpl(new Impl())
{}

OptimumSolverRefiner::OptimumSolverRefiner(const OptimumSolverRefiner& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverRefiner::~OptimumSolverRefiner()
{}

auto OptimumSolverRefiner::operator=(OptimumSolverRefiner other) -> OptimumSolverRefiner&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverRefiner::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverRefiner::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

} // namespace Reaktoro
