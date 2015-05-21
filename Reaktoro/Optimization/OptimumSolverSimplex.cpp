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

#include "OptimumSolverIpfeasible.hpp"

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>
#include <Reaktoro/Optimization/OptimumSolverIpnewton.hpp>

namespace Reaktoro {
namespace {

template<typename Decomposition>
auto solveTranspose(const Decomposition& lu, const Vector& b) -> Vector
{
    const Matrix& LU = lu.matrixLU();

    const auto& L = LU.triangularView<Eigen::UnitLower>();
    const auto& U = LU.triangularView<Eigen::Upper>();
    const auto& P = lu.permutationP();

    Vector x(b.rows());
    x = U.transpose().solve(b);
    x = L.transpose().solve(x);
    x = P.transpose() * x;

    return x;
}

inline auto findFirstNegative(const Vector& vec) -> Index
{
    const Index size = vec.rows();
    Index i = 0; while(i < size and vec[i] >= 0.0) ++i;
    return i;
}

inline auto findMostNegative(const Vector& vec) -> Index
{
    Index i = 0;
    if(vec.rows() == 0) return i;
    else return vec.minCoeff(&i) < 0 ? i : vec.rows();
}

} // namespace

struct OptimumSolverIpfeasible::Impl
{
    OptimumSolverIpnewton ipnewton;

    auto feasible(OptimumProblem problem, OptimumState& state) -> OptimumResult;

    auto simplex(OptimumProblem problem, OptimumState& state) -> OptimumResult;

    auto solve(OptimumProblem problem, OptimumState& state) -> OptimumResult;
};

auto OptimumSolverIpfeasible::Impl::feasible(OptimumProblem problem, OptimumState& state) -> OptimumResult
{
    // The number of variables (n) and equality constraints (m) in the problem
    const unsigned n = problem.c.rows();
    const unsigned m = problem.A.rows();

    const auto& A = problem.A;
    const auto& b = problem.b;
    const auto& lower = problem.lower;
    const auto& upper = problem.upper;

    // Initialise the basic feasible solution of the Phase I problem
    solution.x.resize(n + m);

    auto x1 = solution.x.segment(0, n); // a reference to the first n components of x
    auto x2 = solution.x.segment(n, m); // a reference to the last m components of x

    // Initialise the components of x1 and the partitions L (lower bounds) and U (upper bounds)
    for(unsigned i = 0; i < n; ++i)
    {
        if(std::isfinite(lower[i]))
        {
            x1[i] = lower[i];
            solution.ilower.push_back(i);
        }
        else
        {
            x1[i] = upper[i];
            solution.iupper.push_back(i);
        }
    }

    // Initialise the partition B (basic variables)
    solution.ibasic = range(n, n + m);

    // Initialise the components of x2 (it remains to take the absolute values of x2)
    x2 = b - A * x1;

    // Initialise the Phase I feasibility problem
    LinearProblem feasible_problem(n + m);

    // Initialise the Phase I coefficients `c`
    feasible_problem.c.segment(0, n).setZero();
    feasible_problem.c.segment(n, m).setOnes();

    // Initialise the Phase I equality constraint matrix `A`
    feasible_problem.A.resize(m, n + m);
    feasible_problem.A.leftCols(n) = A;
    for(unsigned i = 0; i < m; ++i)
        feasible_problem.A(i, n + i) = (x2[i] < 0) ? -1.0 : 1.0;

    // Set the components of x2 to its positive values, since `A` has just been assembled
    x2 = x2.array().abs();

    // Initialise the right-hand side vector `b` of the Phase I problem
    feasible_problem.b = b;

    // Initialise the Phase I lower bounds
    feasible_problem.lower.segment(0, n) = lower;   // the lower bounds of the x1 variables
    feasible_problem.lower.segment(n, m).setZero(); // the lower bounds of the x2 variables

    feasible_problem.upper.segment(0, n) = upper;          // the upper bounds of the x1 variables
    feasible_problem.upper.segment(n, m).setConstant(inf); // the upper bounds of the x2 variables

    // Solve the Phase I problem
    simplex(solution, feasible_problem);

    Matrix AB = extractCols(solution.ibasic, feasible_problem.A);

    for(unsigned i = 0; i < m; ++i)
    {
        // The index of the i-th artificial variable
        const Index k = n + i;

        // The local index of the i-th artificial variable in the basic set
        const Index klocal = find(k, solution.ibasic);

        // Check if the i-th artificial variable is in the basic set
        if(klocal < m)
        {
            // Find a variable with index j that can be swapped with the i-th artificial variable
            for(Index j = 0; j < n; ++j)
            {
                // The j-th variable cannot be in the basic set
                if(contains(j, solution.ibasic))
                    continue;

                // Compute the vector `t`
                const Vector t = AB.lu().solve(A.col(j));

                // Check if the i-th component of `t` is non-zero
                if(t[i] != 0.0)
                {
                    // Add the artificial variable `k` to the lower bound set
                    append(k, solution.ilower);  // L' = L + {k}

                    // Add the variable `j` to the basic set and remove the artificial variable `k` from it
                    solution.ibasic[klocal] = j; // B' = B - {k} + {j}

                    // Determine if variable `j` is to be removed from the lower or upper bound set
                    if(solution.x[j] == problem.lower[j])
                        remove(j, solution.ilower); // L' = L - {j}
                    else
                        remove(j, solution.iupper); // U' = U - {j}

                    // Update the basic matrix
                    AB.col(klocal) = feasible_problem.A.col(j);

                    // Since variable j has been found, no need to check other variables
                    break;
                }
            }
        }
    }

    // Remove the artificial variables
    solution.x.conservativeResize(n);
    solution.zl.conservativeResize(n);
    solution.zu.conservativeResize(n);

    auto ilower_end = std::remove_if(solution.ilower.begin(), solution.ilower.end(), [=](Index i) { return i >= n; });
    auto iupper_end = std::remove_if(solution.iupper.begin(), solution.iupper.end(), [=](Index i) { return i >= n; });

    solution.ilower.assign(solution.ilower.begin(), ilower_end);
    solution.iupper.assign(solution.iupper.begin(), iupper_end);
}

auto OptimumSolverIpfeasible::Impl::simplex(OptimumProblem problem, OptimumState& state) -> OptimumResult
{
    const auto& A = problem.A;
    const auto& c = problem.c;
    const auto& lower = problem.lower;
    const auto& upper = problem.upper;

    auto& x  = solution.x;
    auto& y  = solution.y;
    auto& ibasic = solution.ibasic;
    auto& ilower = solution.ilower;
    auto& iupper = solution.iupper;

    const unsigned m = A.rows();
    const unsigned n = A.cols();

    std::sort(ibasic.begin(), ibasic.end());

    Vector xb = extractRows(ibasic, solution.x);

    while(true)
    {
        const unsigned nL = ilower.size();
        const unsigned nU = iupper.size();

        Matrix AB = extractCols(ibasic, A);
        Matrix AL = extractCols(ilower, A);
        Matrix AU = extractCols(iupper, A);

        Vector cB = extractRows(ibasic, c);
        Vector cL = extractRows(ilower, c);
        Vector cU = extractRows(iupper, c);

        Eigen::PartialPivLU<Matrix> lu(AB);

        y = solveTranspose(lu, cB);

        Vector zL =  cL - AL.transpose() * y;
        Vector zU = -cU + AU.transpose() * y;

        const Index qLower = findMostNegative(zL); // the index of the first negative entry in zL
        const Index qUpper = findMostNegative(zU); // the index of the first negative entry in zU

        // Check if all dual variables zL and zU are positive
        if(qLower == nL and qUpper == nU)
        {
            setRows(ibasic, xb, x);
            solution.zl.setZero(n);
            solution.zu.setZero(n);
            setRows(ilower, zL, solution.zl);
            setRows(iupper, zU, solution.zu);
            return;
        }

        // There is a lower bound component with negative dual variable - moving it to the basic set
        if(qLower < nL)
        {
            // The index q of the variable on the lower bound that will be moved to the basic set
            const Index q = ilower[qLower];

            // The local index of q in the set of lower bound variables
            const Index qlocal = qLower;

            // Compute the step vector `t`
            Vector t = lu.solve(A.col(q));

            // Compute the step length `lambda`
            double lambda = upper[q] - lower[q];

            // The code that indicates how the sets B, L, U change
            // 0 : L' = L - {q}       and U' = U + {q}
            // 1 : B' = B + {q} - {p} and L' = L - {q} + {p}
            // 2 : B' = B + {q} - {p} and L' = L - {q} and U' = U + {p}
            int code = 0;

            // The local index p to exit the basic set
            int plocal;

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

            if(not std::isfinite(lambda))
                throw std::runtime_error("***Error***: Unbounded linear programming problem.");

            // The index of the basic variable exiting the basic set
            const Index p = ibasic[plocal];

            switch(code)
            {
            case 0:
                erase(qlocal, ilower);      // L' = L - {q}
                append(q, iupper);          // U' = U + {q}
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
                append(p, iupper);          // U' = U + {p}
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
            Vector t = lu.solve(A.col(q));

            // Compute the step length `lambda`
            double lambda = upper[q] - lower[q];

            // The code that indicates how the sets B, L, U change
            // 0 : U' = U - {q}       and L' = L + {q}
            // 1 : B' = B + {q} - {p} and U' = U - {q} + {p}
            // 2 : B' = B + {q} - {p} and U' = U - {q} and L' = L + {p}
            int code = 0;

            // The local index p to exit the basic set
            int plocal;

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

            if(not std::isfinite(lambda))
                throw std::runtime_error("***Error***: Unbounded linear programming problem.");

            // The index of the basic variable exiting the basic set
            const Index p = ibasic[plocal];

            switch(code)
            {
            case 0:
                erase(qlocal, iupper); // U' = U - {q}
                append(q, ilower);     // L' = L + {q}
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
                append(p, ilower);          // L' = L + {p}
                ibasic[plocal] = q;         // B' = B + {q} - {p}
                xb += lambda * t;           // step to the new vertex
                xb[plocal] = x[q] - lambda; // update the new basic variable
                x[p] = lower[p];            // set the p-th variable to its lower bound
                break;
            }
        }
    }
}

auto OptimumSolverIpfeasible::Impl::solve(OptimumProblem problem, OptimumState& state) -> OptimumResult
{
    feasible(problem, state);
    simplex(problem, state);
}

OptimumSolverIpfeasible::OptimumSolverIpfeasible()
: pimpl(new Impl())
{}

OptimumSolverIpfeasible::OptimumSolverIpfeasible(const OptimumSolverIpfeasible& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverIpfeasible::~OptimumSolverIpfeasible()
{}

auto OptimumSolverIpfeasible::operator=(OptimumSolverIpfeasible other) -> OptimumSolverIpfeasible&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return approximate(problem, state, {});
}

auto OptimumSolverIpfeasible::approximate(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->approximate(problem, state, options);
}

} // namespace Reaktoro
