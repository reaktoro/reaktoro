// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "OptimumSolverSimplex.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/Utils.hpp>

namespace Reaktoro {
namespace {

template<typename Decomposition>
auto solveTranspose(const Decomposition& lu, VectorXdConstRef b) -> VectorXd
{
    const MatrixXd& LU = lu.matrixLU();

    const auto& L = LU.triangularView<Eigen::UnitLower>();
    const auto& U = LU.triangularView<Eigen::Upper>();
    const auto& P = lu.permutationP();

    VectorXd x(b.rows());
    x = U.transpose().solve(b);
    x = L.transpose().solve(x);
    x = P.transpose() * x;

    return x;
}

inline auto findFirstNegative(VectorXdConstRef vec) -> Index
{
    const Index size = vec.rows();
    Index i = 0; while(i < size && vec[i] >= 0.0) ++i;
    return i;
}

inline auto findMostNegative(VectorXdConstRef vec) -> Index
{
    Index i = 0;
    if(vec.rows() == 0) return i;
    else return vec.minCoeff(&i) < 0 ? i : vec.rows();
}

template<typename T>
inline void erase(const Index& i, std::vector<T>& vec)
{
    vec.erase(vec.begin() + i);
}

template<typename T, typename U>
inline void remove(const U& element, std::vector<T>& vec)
{
    std::remove(vec.begin(), vec.end(), element);
}

} // namespace

struct OptimumSolverSimplex::Impl
{
    Indices ibasic, ilower, iupper;

    auto feasible(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    auto simplex(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverSimplex::Impl::feasible(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    OptimumResult result;

    // The number of variables (n) and equality constraints (m) in the problem
    const Index n = problem.c.rows();
    const Index m = problem.A.rows();

    // Define some auxiliary references
    const auto& A = problem.A;
    const auto& b = problem.b;
    const auto& lower = problem.l;
    const auto& upper = problem.u;

    // Initialise the basic feasible solution of the Phase I problem
    state.x.resize(n + m);

    auto x1 = state.x.segment(0, n); // a reference to the first n components of x
    auto x2 = state.x.segment(n, m); // a reference to the last m components of x

    // Initialise the components of x1 and the partitions L (lower bounds) and U (upper bounds)
    ibasic.clear();
    ilower.clear();
    iupper.clear();
    for(unsigned i = 0; i < n; ++i)
    {
        if(std::isfinite(lower[i]))
        {
            x1[i] = lower[i];
            ilower.push_back(i);
        }
        else
        {
            x1[i] = upper[i];
            iupper.push_back(i);
        }
    }

    // Initialise the partition B (basic variables)
    ibasic = range(n, n + m);

    // Initialise the components of x2 (it remains to take the absolute values of x2)
    x2 = b - A * x1;

    // Initialise the Phase I feasibility problem
    OptimumProblem feasible_problem;

    // Initialise the Phase I coefficients `c`
    feasible_problem.c.resize(n + m);
    feasible_problem.c.segment(0, n).setZero();
    feasible_problem.c.segment(n, m).setOnes();

    // Initialise the Phase I equality constraint matrix `A`
    feasible_problem.A = zeros(m, n + m);
    feasible_problem.A.leftCols(n) = A;
    for(unsigned i = 0; i < m; ++i)
        feasible_problem.A(i, n + i) = (x2[i] < 0) ? -1.0 : 1.0;

    // Set the components of x2 to its positive values, since `A` has just been assembled
    x2 = x2.array().abs();

    // Initialise the right-hand side vector `b` of the Phase I problem
    feasible_problem.b = b;

    // Initialise the Phase I lower bounds
    feasible_problem.l.resize(n + m);
    feasible_problem.l.segment(0, n) = lower;   // the lower bounds of the x1 variables
    feasible_problem.l.segment(n, m).setZero(); // the lower bounds of the x2 variables

    feasible_problem.u.resize(n + m);
    feasible_problem.u.segment(0, n) = upper;                 // the upper bounds of the x1 variables
    feasible_problem.u.segment(n, m).setConstant(infinity()); // the upper bounds of the x2 variables

    // Solve the Phase I problem
    result = simplex(feasible_problem, state, options);

    // Check if a basic feasible solution exists by checking the sum of artificial variables
    Assert(x2.sum() <= 1e-16 * x1.sum() && result.succeeded,
        "Failed to calculate a basic feasible solution.",
        "The provided constraints result in an infeasible problem.");

    // Check if the basic feasible solution has an artificial variable
    for(unsigned k = n; k < n + m; ++k)
        if(contained(k, ibasic))
            RuntimeError("A basic feasible solution was calculated, but it contains artificial variables.",
                "The provided constraints might be linearly dependent.");

    // Remove the artificial variables
    state.x.conservativeResize(n);
    state.z.conservativeResize(n);
    state.w.conservativeResize(n);

    auto ilower_end = std::remove_if(ilower.begin(), ilower.end(), [=](Index i) { return i >= n; });
    auto iupper_end = std::remove_if(iupper.begin(), iupper.end(), [=](Index i) { return i >= n; });

    ilower.assign(ilower.begin(), ilower_end);
    iupper.assign(iupper.begin(), iupper_end);

    return result;
}

auto OptimumSolverSimplex::Impl::simplex(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Start timing the calculation
    Time begin = time();

    // The result of the calculation
    OptimumResult result;

    // Define some auxiliary references to problem variables
    const auto& A = problem.A;
    const auto& c = problem.c;
    const auto& lower = problem.l;
    const auto& upper = problem.u;
    const unsigned m = A.rows();
    const unsigned n = A.cols();

    // Define some auxiliary references to state variables
    auto& x = state.x;
    auto& y = state.y;
    auto& z = state.z;
    auto& w = state.w;
    auto& f = state.f;

    // Define some auxiliary references to result variables
    auto& iterations = result.iterations;

    // Define auxiliary references to general options
    const auto maxiters = options.max_iterations;

    std::sort(ibasic.begin(), ibasic.end());

    VectorXd xb = rows(state.x, ibasic);

    for(iterations = 1; iterations <= maxiters; ++iterations)
    {
        const unsigned nL = ilower.size();
        const unsigned nU = iupper.size();

        MatrixXd AB = cols(A, ibasic);
        MatrixXd AL = cols(A, ilower);
        MatrixXd AU = cols(A, iupper);

        VectorXd cB = rows(c, ibasic);
        VectorXd cL = rows(c, ilower);
        VectorXd cU = rows(c, iupper);

        Eigen::PartialPivLU<MatrixXd> lu(AB);

        y = solveTranspose(lu, cB);

        VectorXd zL =  cL - AL.transpose() * y;
        VectorXd zU = -cU + AU.transpose() * y;

        const Index qLower = findMostNegative(zL); // the index of the first negative entry in zL
        const Index qUpper = findMostNegative(zU); // the index of the first negative entry in zU

        // Check if all dual variables zL and zU are positive
        if(qLower == nL && qUpper == nU)
        {
            rows(x, ibasic) = xb;
            z.setZero(n);
            w.setZero(n);
            rows(z, ilower) = zL;
            rows(w, iupper) = zU;
            break;
        }

        // There is a lower bound component with negative dual variable - moving it to the basic set
        if(qLower < nL)
        {
            // The index q of the variable on the lower bound that will be moved to the basic set
            const Index q = ilower[qLower];

            // The local index of q in the set of lower bound variables
            const Index qlocal = qLower;

            // Compute the step vector `t`
            VectorXd t = lu.solve(A.col(q));

            // Compute the step length `lambda`
            double lambda = upper[q] - lower[q];

            // The code that indicates how the sets B, L, U change
            // 0 : L' = L - {q}       and U' = U + {q}
            // 1 : B' = B + {q} - {p} and L' = L - {q} + {p}
            // 2 : B' = B + {q} - {p} and L' = L - {q} and U' = U + {p}
            int code = 0;

            // The local index p to exit the basic set
            int plocal = -1;

            for(unsigned i = 0; i < m; ++i)
            {
                if(t[i] <= 0.0) continue; // skip for negative components of `t`

                // Compute the step to bind the i-th basic variable to its lower bound
                const double trial = (xb[i] - lower[ibasic[i]])/t[i];

                if(trial < lambda)
                {
                    lambda = trial;
                    plocal = i;
                    code   = 1;
                }
            }

            for(unsigned i = 0; i < m; ++i)
            {
                if(t[i] >= 0.0) continue; // skip for positive components of `t`

                // Compute the step to bind the i-th basic variable to its upper bound
                const double trial = (xb[i] - upper[ibasic[i]])/t[i];

                if(trial < lambda)
                {
                    lambda = trial;
                    plocal = i;
                    code   = 2;
                }
            }

            // Check if the linear programming programming is bounded
            Assert(std::isfinite(lambda),
                "Cannot proceed with the simplex minimization.",
                "The linear programming problem is unbounded.");

            // The index of the basic variable exiting the basic set
            const Index p = ibasic[plocal];

            switch(code)
            {
            case 0:
                erase(qlocal, ilower);      // L' = L - {q}
                iupper.push_back(q);        // U' = U + {q}
                x[q] = upper[q];            // set the q-th variable to its upper bound
                break;

            case 1:
                ilower[qlocal] = p;         // L' = L - {q} + {p}
                ibasic[plocal] = q;         // B' = B + {q} - {p}
                xb -= lambda * t;           // step to the new vertex
                xb[plocal] = x[q] + lambda; // update the new basic variable
                x[p] = lower[p];            // set the p-th variable to its lower bound
                break;

            case 2:
                erase(qlocal, ilower);      // L' = L - {q}
                iupper.push_back(p);        // U' = U + {p}
                ibasic[plocal] = q;         // B' = B + {q} - {p}
                xb -= lambda * t;           // step to the new vertex
                xb[plocal] = x[q] + lambda; // update the new basic variable
                x[p] = upper[p];            // set the p-th variable to its upper bound
                break;
            }
        }

        // There is an upper bound component with negative dual variable - moving it to the basic set
        else
        {
            // The index q of the variable on the upper bound that will be moved to the basic set
            const Index q = iupper[qUpper];

            // The local index of q in the set of upper bound variables
            const Index qlocal = qUpper;

            // Compute the step vector `t`
            VectorXd t = lu.solve(A.col(q));

            // Compute the step length `lambda`
            double lambda = upper[q] - lower[q];

            // The code that indicates how the sets B, L, U change
            // 0 : U' = U - {q}       and L' = L + {q}
            // 1 : B' = B + {q} - {p} and U' = U - {q} + {p}
            // 2 : B' = B + {q} - {p} and U' = U - {q} and L' = L + {p}
            int code = 0;

            // The local index p to exit the basic set
            int plocal = -1;

            for(unsigned i = 0; i < m; ++i)
            {
                if(t[i] <= 0.0) continue; // skip for negative components of `t`

                // Compute the step to bind the i-th basic variable to its upper bound
                const double trial = (upper[ibasic[i]] - xb[i])/t[i];

                if(trial < lambda)
                {
                    lambda = trial;
                    plocal = i;
                    code   = 1;
                }
            }

            for(unsigned i = 0; i < m; ++i)
            {
                if(t[i] >= 0.0) continue; // skip for positive components of `t`

                // Compute the step to bind the i-th basic variable to its lower bound
                const double trial = (lower[ibasic[i]] - xb[i])/t[i];

                if(trial < lambda)
                {
                    lambda = trial;
                    plocal = i;
                    code   = 2;
                }
            }

            // Check if the linear programming programming is bounded
            Assert(std::isfinite(lambda),
                "Cannot proceed with the simplex minimization.",
                "The linear programming problem is unbounded.");

            // The index of the basic variable exiting the basic set
            const Index p = ibasic[plocal];

            switch(code)
            {
            case 0:
                erase(qlocal, iupper); // U' = U - {q}
                ilower.push_back(q);   // L' = L + {q}
                x[q] = lower[q];       // set the q-th variable to its lower bound
                break;

            case 1:
                iupper[qlocal] = p;         // U' = U - {q} + {p}
                ibasic[plocal] = q;         // B' = B + {q} - {p}
                xb += lambda * t;           // step to the new vertex
                xb[plocal] = x[q] - lambda; // update the new basic variable
                x[p] = upper[p];            // set the p-th variable to its upper bound
                break;

            case 2:
                erase(qlocal, iupper);      // U' = U - {q}
                ilower.push_back(p);        // L' = L + {p}
                ibasic[plocal] = q;         // B' = B + {q} - {p}
                xb += lambda * t;           // step to the new vertex
                xb[plocal] = x[q] - lambda; // update the new basic variable
                x[p] = lower[p];            // set the p-th variable to its lower bound
                break;
            }
        }
    }

    // Set the state of the objective result
    f.val = dot(c, x);
    f.grad = c;

    // Check if calculation was successful
    result.succeeded = iterations < maxiters;

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

auto OptimumSolverSimplex::Impl::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    // Solve the regularized optimization problem
    auto result = feasible(problem, state, options);
    result += simplex(problem, state, options);
    return result;
}

OptimumSolverSimplex::OptimumSolverSimplex()
: pimpl(new Impl())
{}

OptimumSolverSimplex::OptimumSolverSimplex(const OptimumSolverSimplex& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverSimplex::~OptimumSolverSimplex()
{}

auto OptimumSolverSimplex::operator=(OptimumSolverSimplex other) -> OptimumSolverSimplex&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverSimplex::feasible(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->feasible(problem, state, options);
}

auto OptimumSolverSimplex::simplex(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->simplex(problem, state, options);
}

auto OptimumSolverSimplex::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    OptimumOptions options;
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverSimplex::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverSimplex::dxdp(VectorXdConstRef dgdp, VectorXdConstRef dbdp) -> VectorXd
{
    RuntimeError("Could not calculate the sensitivity of the optimal solution with respect to parameters.",
        "The method OptimumSolverSimplex::dxdp has not been implemented yet.");
    return {};
}

auto OptimumSolverSimplex::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverSimplex(*this);
}

} // namespace Reaktoro
