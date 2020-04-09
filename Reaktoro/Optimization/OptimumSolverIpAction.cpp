// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "OptimumSolverIpAction.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Outputter.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Math/LU.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/KktSolver.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {
namespace {

struct LinearSystemSolverDiagonal
{
    VectorXd A1, A2, invA1, a1, a2, q, u;
    MatrixXd B1, B2, C1, C2, Q;
    Indices ipivot, inonpivot;
    Index n, m;
    Eigen::PartialPivLU<MatrixXd> lu;

    /// Decompose the matrix [A B ; C I]
    auto decompose(VectorXdConstRef A, MatrixXdConstRef B, MatrixXdConstRef C) -> void
    {
        n = A.rows();
        m = B.cols();
        ipivot.clear();
        inonpivot.clear();
        ipivot.reserve(n);
        inonpivot.reserve(n);
        for(unsigned i = 0; i < n; ++i)
            if(A[i] > norminf(C.col(i))) ipivot.push_back(i);
            else inonpivot.push_back(i);

        const unsigned n2 = inonpivot.size();

        A1 = rows(A, ipivot);
        A2 = rows(A, inonpivot);
        invA1 = 1/A1;
        B1 = rows(B, ipivot);
        B2 = rows(B, inonpivot);
        C1 = cols(C, ipivot);
        C2 = cols(C, inonpivot);

        Q = zeros(n2 + m, n2 + m);
        Q.topLeftCorner(n2, n2).diagonal() = A2;
        Q.topRightCorner(n2, m) = B2;
        Q.bottomLeftCorner(m, n2) = C2;
        Q.bottomRightCorner(m, m) = identity(m, m);
        Q.bottomRightCorner(m, m) -= C1*diag(invA1)*B1;

        lu.compute(Q);
    }

    /// Solve the linear system
    auto solve(VectorXdConstRef a, VectorXdConstRef b, VectorXdRef x, VectorXdRef y) -> bool
    {
        const unsigned n2 = inonpivot.size();

        a1 = rows(a, ipivot);
        a2 = rows(a, inonpivot);

        q.resize(n2 + m);
        rows(q, 0, n2) = a2;
        rows(q, n2, m) = b - C1*diag(invA1)*a1;

        u = lu.solve(q);

        if(!u.allFinite())
            u = Q.fullPivLu().solve(q);

        y = rows(u, n2, m);

        x.resize(n);
        rows(x, ipivot) = invA1 % (a1 - B1*y);
        rows(x, inonpivot) = rows(u, 0, n2);

        return x.allFinite() && y.allFinite();
    }
};

} // namespace

struct OptimumSolverIpAction::Impl
{
    /// The right-hand side of the KKT equations
    LinearSystemSolverDiagonal lssd;

    /// The outputter instance
    Outputter outputter;

    /// The weighted LU decomposition of the coefficient matrix A
    LU lu;

    /// The coefficient matrix of the secondary variables
    MatrixXd S;

    /// The kernel matrix
    MatrixXd K;

    /// The indices of the primary variables
    Indices iP;

    /// The indices of the secondary variables
    Indices iS;

    /// The trial iterate x
    VectorXd xtrial;

    /// The Newton steps for variables x and z
    VectorXd dx, dz;

    /// The Newton steps for primary and secondary variables
    VectorXd dxP, dxS;

    /// The residual vectors in the Newton step equation
    VectorXd r1, r2;

    /// The diagonal Hessian of the Lagrange function and its primary and secondary components
    VectorXd H, HP, HS;

    /// The auxiliary matrix StHP = tr(S)*HP
    MatrixXd StHP;

    /// The names of the secondary variables
    std::vector<std::string> snames;

    /// Solve the optimization problem.
    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
    {
        // Start timing the calculation
        Time begin = time();

        // The result of the calculation
        OptimumResult result;

        // Finish the calculation if the problem has no variable
        if(problem.n == 0)
        {
            state = OptimumState();
            result.succeeded = true;
            result.time = elapsed(begin);
            return result;
        }

        // Initialize the outputter instance
        outputter = Outputter();
        outputter.setOptions(options.output);

        // Define some auxiliary references to variables
        auto& x = state.x;
        auto& y = state.y;
        auto& z = state.z;
        auto& f = state.f;

        // The number of variables and equality constraints
        const auto& n = problem.A.cols();
        const auto& m = problem.A.rows();

        // The components of the equality constraints
        MatrixXd A = problem.A;
        VectorXd b = problem.b;

        // Define auxiliary references to general options
        const auto tol = options.tolerance;
        const auto maxiters = options.max_iterations;

        // Define some auxiliary references to IpAction parameters
        const auto mu = options.ipaction.mu;
        const auto tau = options.ipaction.tau;

        // Define some auxiliary references to result variables
        auto& error = result.error;
        auto& iterations = result.iterations;
        auto& succeeded = result.succeeded = false;

        // The regularization parameters delta and gamma
        auto gamma = options.regularization.gamma;
        auto delta = options.regularization.delta;

        // Ensure the initial guesses for `x` and `y` have adequate dimensions
        if(x.size() != n) x = zeros(n);
        if(y.size() != m) y = zeros(m);
        if(z.size() != n) z = zeros(n);

        // Ensure the initial guesses for `x` and `z` are inside the feasible domain
        x = (x.array() > 0.0).select(x, mu);
        z = (z.array() > 0.0).select(z, 1.0);

        // The optimality, feasibility, centrality and total error variables
        double errorf, errorh, errorc;

        // The function that outputs the header and initial state of the solution
        auto output_initial_state = [&]()
        {
            if(!options.output.active) return;

            // Initialize the names of the secondary variables
            snames = extract(options.output.xnames, iS);

            outputter.addEntry("Iteration");
            outputter.addEntries(options.output.xprefix, n, options.output.xnames);
            outputter.addEntries(options.output.yprefix, m, options.output.ynames);
            outputter.addEntries(options.output.zprefix, n, options.output.znames);
            outputter.addEntries("r", snames.size(), snames);
            outputter.addEntry("f(x)");
            outputter.addEntry("Error");
            outputter.addEntry("Optimality");
            outputter.addEntry("Feasibility");
            outputter.addEntry("Centrality");

            outputter.outputHeader();
            outputter.addValue(iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValues(abs(r1));
            outputter.addValue(f.val);
            outputter.addValue(error);
            outputter.addValue(errorf);
            outputter.addValue(errorh);
            outputter.addValue(errorc);
            outputter.outputState();
        };

        // The function that outputs the current state of the solution
        auto output_state = [&]()
        {
            if(!options.output.active) return;

            outputter.addValue(iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(z);
            outputter.addValues(abs(r1));
            outputter.addValue(f.val);
            outputter.addValue(error);
            outputter.addValue(errorf);
            outputter.addValue(errorh);
            outputter.addValue(errorc);
            outputter.outputState();
        };

        auto initialize_decomposition = [&]()
        {
            // Initialize the weighted scaling vector `W`
            VectorXd W = abs(x);
            const double threshold = 1e-10 * (max(W) + 1);
            W = (W.array() > threshold).select(W, threshold);

            // Compute the weighted LU decomposition of the coefficient matrix A
            lu.compute(A, W);

            // Auxiliary references to LU factors
            const auto& r = lu.rank;
            const auto& P = lu.P;
            const auto& Q = lu.Q;
            const auto& U = lu.U;
            const auto& L = lu.L.topLeftCorner(r, r).triangularView<Eigen::Lower>();

            // Initialize the indices of the primary variables
            iP = Indices(Q.indices().data(), Q.indices().data() + m);

            // Initialize the indices of the secondary variables
            iS = Indices(Q.indices().data() + m, Q.indices().data() + n);

            // References to the U1 and U2 part of U = [U1 U2]
            const auto& U1 = U.topLeftCorner(r, r).triangularView<Eigen::Upper>();
            const auto& U2 = U.topRightCorner(r, n - r);

            // Initialize the coefficient matrix S
            S = U2;
            S = U1.solve(S);

            // Clean the matrix S from round-off errors if composed by rational numbers
            if(options.regularization.max_denominator)
                cleanRationalNumbers(S, options.regularization.max_denominator);

            // Initialize the kernel matrix K
            K.resize(n - r, n);
            K.leftCols(r) = tr(S);
            K.rightCols(n - r) = -identity(n - r, n - r);
            K = K * Q.inverse();

            // Clean the kernel matrix from round-off errors
            if(options.regularization.max_denominator)
                cleanRationalNumbers(K, options.regularization.max_denominator);

            // Initialize the transformed equality constraints matrix A
            A.resize(r, n);
            A.leftCols(r) = identity(r, r);
            A.rightCols(n - r) = S;
            A = A * Q.inverse();

            // Clean the transformed equality constraints matrix A from round-off errors
            if(options.regularization.max_denominator)
                cleanRationalNumbers(A, options.regularization.max_denominator);

            // Initialize the transformed equality constraints vector b
            b = P * b;
            b.conservativeResize(r);
            b = L.solve(b);
            b = U1.solve(b);
        };

        // Return true if the result of a calculation failed
        auto failed = [&](bool succeeded)
        {
            return !succeeded;
        };

        // The function that computes the current error norms
        auto update_residuals = [&]()
        {
            // Compute the right-hand side vectors of the KKT equation
            r1 = -K*(f.grad - mu/x);
            r2 = -(A*x + delta*delta*y - b);

            // Calculate the optimality, feasibility and centrality errors
            errorf = norminf(r1);
            errorh = norminf(r2);
            errorc = norminf(x % z - mu);
            error = std::max({errorf, errorh, errorc});
        };

        // The function that initialize the state of some variables
        auto initialize = [&]()
        {
            // Initialize xtrial
            xtrial.resize(n);

            // Evaluate the objective function
            f = problem.objective(x);

            // Update the residuals of the calculation
            update_residuals();
        };

        // The function that computes the Action step
        auto compute_newton_step_diagonal = [&]()
        {
            // Calculate the diagonal Hessian of the Lagrange function
            H = f.hessian.diagonal + z/x + gamma*gamma*ones(n);

            // Extract the primary and secondary components
            HP = rows(H, iP);
            HS = rows(H, iS); HS *= -1;

            // Compute the auxiliary matrix StHP = tr(S)*HP
            StHP = tr(S)*diag(HP);

            // Update the decomposition of the linear system
            lssd.decompose(HS, StHP, S);

            // Compute the steps for the primary and secondary variables
            bool successful = lssd.solve(r1, r2, dxS, dxP);

            // Assemble the Newton step dx from the primary and secondary steps dxP and dxS
            dx.resize(n);
            rows(dx, iP) = dxP;
            rows(dx, iS) = dxS;

            // Compute the dz step using dx
            dz = (mu - z % x - z % dx)/x;

            // Return true if he calculation succeeded
            return successful;
        };

        // The aggressive mode for updating the iterates
        auto update_iterates_aggressive = [&]()
        {
            // Calculate the current trial iterate for x
            for(int i = 0; i < n; ++i)
                xtrial[i] = (x[i] + dx[i] > 0.0) ?
                    x[i] + dx[i] : x[i]*(1.0 - tau);

            // Evaluate the objective function at the trial iterate
            f = problem.objective(xtrial);

            // Initialize the step length factor
            double alpha = fractionToTheBoundary(x, dx, tau);

            // The number of tentatives to find a trial iterate that results in finite objective result
            unsigned tentatives = 0;

            // Repeat until f(xtrial) is finite
            while(!isfinite(f) && ++tentatives < 10)
            {
                // Calculate a new trial iterate using a smaller step length
                xtrial = x + alpha * dx;

                // Evaluate the objective function at the trial iterate
                f = problem.objective(xtrial);

                // Decrease the current step length
                alpha *= 0.5;
            }

            // Return false if xtrial could not be found s.t. f(xtrial) is finite
            if(tentatives == 10)
                return false;

            // Update the iterate x from xtrial
            x = xtrial;

            // Update the z-Lagrange multipliers
            for(int i = 0; i < n; ++i)
                z[i] += (z[i] + dz[i] > 0.0) ?
                    dz[i] : -tau * z[i];

            // Compute the Lagrange multipliers y only if regularization param delta is non-zero
            if(delta) y = lu.trsolve(f.grad - z + gamma*gamma*x);

            // Return true as found xtrial results in finite f(xtrial)
            return true;
        };

        // The conservative mode for updating the iterates
        auto update_iterates_convervative = [&]()
        {
            // Initialize the step length factor
            double alphax = fractionToTheBoundary(x, dx, tau);
            double alphaz = fractionToTheBoundary(z, dz, tau);
            double alpha = alphax;

            // The number of tentatives to find a trial iterate that results in finite objective result
            unsigned tentatives = 0;

            // Repeat until a suitable xtrial iterate if found such that f(xtrial) is finite
            for(; tentatives < 10; ++tentatives)
            {
                // Calculate the current trial iterate for x
                xtrial = x + alpha * dx;

                // Evaluate the objective function at the trial iterate
                f = problem.objective(xtrial);

                // Leave the loop if f(xtrial) is finite
                if(isfinite(f))
                    break;

                // Decrease alpha in a hope that a shorter step results f(xtrial) finite
                alpha *= 0.01;
            }

            // Return false if xtrial could not be found s.t. f(xtrial) is finite
            if(tentatives == 10)
                return false;

            // Update the iterate x from xtrial
            x = xtrial;

            // Update the z-Lagrange multipliers
            z += alphaz * dz;

            // Compute the Lagrange multipliers y only if regularization param delta is non-zero
            if(delta) y = lu.trsolve(f.grad - z + gamma*gamma*x);

            // Return true as found xtrial results in finite f(xtrial)
            return true;
        };

        // The function that performs an update in the iterates
        auto update_iterates = [&]()
        {
            switch(options.ipnewton.step)
            {
            case Aggressive: return update_iterates_aggressive();
            default: return update_iterates_convervative();
            }
        };

        auto converged = [&]()
        {
            return error < tol;
        };

        initialize_decomposition();
        initialize();
        output_initial_state();

        for(iterations = 1; iterations <= maxiters && !succeeded; ++iterations)
        {
            if(failed(compute_newton_step_diagonal()))
                break;
            if(failed(update_iterates()))
                break;
            if((succeeded = converged()))
                break;
            update_residuals();
            output_state();
        }

        // Output a final header
        outputter.outputHeader();

        // Finish timing the calculation
        result.time = elapsed(begin);

        return result;
    }

    /// Calculate the sensitivity of the optimal solution with respect to parameters.
    auto dxdp(VectorXdConstRef dgdp, VectorXdConstRef dbdp) -> MatrixXd
    {
        // Initialize the right-hand side of the KKT equations
        r1.noalias() = -dgdp;
        r2.noalias() =  dbdp;

        // Solve the linear system equations to get the sensitivities
        lssd.solve(r1, r2, dxS, dxP);

        // Transfer primary and secondary dxP and dxS to dx
        rows(dx, iP) = dxP;
        rows(dx, iS) = dxS;

        // Return the calculated sensitivity vector
        return dx;
    }
};

OptimumSolverIpAction::OptimumSolverIpAction()
: pimpl(new Impl())
{}

OptimumSolverIpAction::OptimumSolverIpAction(const OptimumSolverIpAction& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpAction::~OptimumSolverIpAction()
{}

auto OptimumSolverIpAction::operator=(OptimumSolverIpAction other) -> OptimumSolverIpAction&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpAction::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverIpAction::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverIpAction::dxdp(VectorXdConstRef dgdp, VectorXdConstRef dbdp) -> VectorXd
{
    return pimpl->dxdp(dgdp, dbdp);
}

auto OptimumSolverIpAction::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverIpAction(*this);
}

} // namespace Reaktoro
