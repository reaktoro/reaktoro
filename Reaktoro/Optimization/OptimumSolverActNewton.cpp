// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "OptimumSolverActNewton.hpp"

// Eigen includes
#include <Reaktoro/Math/Eigen/Cholesky>
#include <Reaktoro/Math/Eigen/LU>

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
namespace {

//auto multiKahanSum(const Matrix& A, const Vector& x, Vector& res) -> void
//{
//    res = zeros(A.rows());
//    for(int i = 0; i < A.rows(); ++i)
//    {
//        double c = 0.0;
//        for(int j = 0; j < A.cols(); ++j)
//        {
//            volatile double y = A(i, j)*x[j] - c;
//            volatile double t = res[i] + y;
//            c = (t - res[i]) - y;
//            res[i] = t;
//        }
//    }
//}

auto stepLengthToBound(const Vector& p, const Vector& dp, Index& ilimiting) -> double
{
    ilimiting = p.size();
    double alpha_max = 1.0;
    for(int i = 0; i < p.size(); ++i)
    {
        if(dp[i] < 0.0 && p[i] != 0.0)
        {
            const double alpha_trial = -p[i]/dp[i];
            if(alpha_trial < alpha_max)
            {
                alpha_max = alpha_trial;
                ilimiting = i;
            }
        }
    }

    return alpha_max;
}

auto erase(Vector& x, Index i) -> void
{
    std::copy(x.data() + i + 1, x.data() + x.size(), x.data() + i);
    x.conservativeResize(x.size() - 1);
}

auto erase(Indices& x, Index i) -> void
{
    x.erase(x.begin() + i);
}

} // namespace

struct OptimumSolverActNewton::Impl
{
    KktVector rhs;
    KktSolution sol;
    KktSolver kkt;

    Indices L, F, E;
    Vector gF, gL, gE;
    Vector xF, zF, xE, zL;
    Matrix AF, AL, AE;
    Hessian HF;

    Eigen::LLT<Matrix> llt;

    Outputter outputter;

    auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;
};

auto OptimumSolverActNewton::Impl::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
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
    auto& f = state.f;

    const double delta = options.regularization.delta;
    const double gamma = options.regularization.gamma;

    // The number of variables and equality constraints
    const auto n = problem.A.cols();
    const auto m = problem.A.rows();

    const auto& A = problem.A;
    const auto& b = problem.b;
    const auto& l = problem.l;

    Vector h;

    const double threshold = options.actnewton.threshold;

    // Ensure `x` has dimension `n` and its components do not violate the bounds
    x.resize(n);
    x.noalias() = (x.array() > l.array()).select(x, l);

    // Ensure `y` has dimension `m` and proper initial values
    if(y.size() != m) y = zeros(m);

    // Initialize the set of free variables
    for(int i = 0; i < n; ++i)
        if(x[i] == l[i])
            L.push_back(i);
        else if(x[i] > threshold)
            F.push_back(i);
        else
            E.push_back(i);

    // Initialize the submatrices AF and AL from A
    AF = cols(A, F);
    AL = cols(A, L);
    AE = cols(A, E);
    xF = rows(x, F);
    xE = rows(x, E);

    // The alpha step sizes used to restric the steps inside the feasible domain
    double alphax, alphaz, alpha;

    // The optimality, feasibility, centrality and total error variables
    double errorf, errorh, error;

    Index iglimiting;

    // The function that outputs the header and initial state of the solution
    auto output_header = [&]()
    {
        if(!options.output.active) return;

        outputter.addEntry("iter");
        outputter.addEntries(options.output.xprefix, n, options.output.xnames);
        outputter.addEntries(options.output.yprefix, m, options.output.ynames);
        outputter.addEntries(options.output.zprefix, n, options.output.znames);
        outputter.addEntries("r", n, options.output.xnames);
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("errorf");
        outputter.addEntry("errorh");
        outputter.addEntry("error");
        outputter.addEntry("alpha");
        outputter.addEntry("alphax");
        outputter.addEntry("alphaz");

        Vector r = abs(f.grad - tr(A)*y);

        outputter.outputHeader();
        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValues(r);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue("---");
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

        Vector r = abs(f.grad - tr(A)*y);

        outputter.addValue(result.iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(z);
        outputter.addValues(r);
        outputter.addValue(f.val);
        outputter.addValue(norminf(h));
        outputter.addValue(errorf);
        outputter.addValue(errorh);
        outputter.addValue(error);
        outputter.addValue(alpha);
        outputter.addValue(alphax);
        outputter.addValue(alphaz);
        outputter.outputState();
    };

    // The function that updates the objective and constraint state
    auto update_state = [&]()
    {
        rows(x, E) = xE;
        rows(x, F) = xF;
        rows(x, L) = rows(l, L);

        f = problem.objective(x);
        h = A*x - b;

        if(y.norm() == 0.0)
        {
            gF = rows(f.grad, F);
            y = tr(AF).fullPivLu().solve(gF);
        }

        // Update the gradient subvector corresponding to the interior variables
        gL = rows(f.grad, L);
        zL = gL - tr(AL)*y;

        rows(z, L) = zL;

        if(L.size())
        {
            Index iminz;
            const double minz = zL.minCoeff(&iminz);

            if(minz < 0)
            {
                outputter.outputMessage("Adding ", options.output.xnames[L[iminz]], " to the set of free variables...\n");
                F.push_back(L[iminz]);
                erase(L, iminz);
                xF.conservativeResize(F.size());
                xF[F.size()-1] = l[F.back()];
                AF = cols(A, F);
                AL = cols(A, L);
            }
        }

        for(Index i = 0; i < E.size(); ++i)
        {
            if(xE[i] > threshold)
            {
                outputter.outputMessage("Adding ", options.output.xnames[E[i]], " to the set of free variables...\n");
                F.push_back(E[i]);
                xF.conservativeResize(F.size());
                xF[F.size()-1] = xE[i];
                erase(E, i);
                erase(xE, i);
            }
        }

        AE = cols(A, E);
        AF = cols(A, F);

        gE = rows(f.grad, E);
        gF = rows(f.grad, F);

        // Update the Hessian submatrix corresponding to the interior variables
        HF.mode = f.hessian.mode;
        switch(f.hessian.mode)
        {
        case Hessian::Dense:
            HF.dense = submatrix(f.hessian.dense, F, F);
            break;
        case Hessian::Diagonal:
            HF.diagonal = rows(f.hessian.diagonal, F); break;
        default:
            RuntimeError("Could not solve the optimization problem with given Hessian.",
                "OptimumSolverActNewton only accepts `Dense` or `Diagonal` Hessian matrices.");
        }
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
    auto compute_newton_step = [&]()
    {
        zF = zeros(F.size());
        KktMatrix lhs(HF, AF, xF, zF);

        Vector CF = delta*delta*ones(F.size());
        kkt.decompose(lhs);

        // Compute the right-hand side vectors of the KKT equation
        rhs.rx.noalias() = -(gF - tr(AF)*y + gamma*gamma*xF);
        rhs.ry.noalias() = -(h + delta*delta*y);
        rhs.rz.noalias() = zeros(F.size());

        // Compute `dx` and `dy` by solving the KKT equation
        kkt.solve(rhs, sol);

        // Update the time spent in linear systems
        result.time_linear_systems += kkt.result().time_solve;
        result.time_linear_systems += kkt.result().time_decompose;
    };

    // Return true if the function `compute_newton_step` failed
    auto compute_newton_step_failed = [&]()
    {
        const bool dx_finite = sol.dx.allFinite();
        const bool dy_finite = sol.dy.allFinite();
        const bool dz_finite = sol.dz.allFinite();
        const bool all_finite = dx_finite && dy_finite && dz_finite;
        return !all_finite;
    };

    // The function that performs an update in the iterates
    auto update_iterates = [&]()
    {
        Index ilimiting;
        Vector lE = rows(l, E);
        Vector lF = rows(l, F);
        alpha = stepLengthToBound(xF-lF, sol.dx, ilimiting);

        iglimiting = F[ilimiting];

        xF += alpha * sol.dx;
        y += sol.dy;

        xF = (xF.array() > lF.array()).select(xF, lF);
        xE = (xE.array() > lE.array()).select(xE, lE);

        rows(x, F) = xF;
        rows(x, E) = xE;

        // Check if there is a limiting variable that should become active on the bound
        double zlimiting;
        std::string name;

        if(ilimiting < F.size())
        {
            zlimiting = gF[ilimiting] - dot(AF.col(ilimiting), y);
            name = options.output.xnames[F[ilimiting]];
        }

        if(ilimiting < F.size() && zlimiting > 0.0)
        {
            outputter.outputMessage("Adding ", name, " to the set of blocked variables...\n");
            Index iF = F[ilimiting];
            L.push_back(iF);
            erase(F, ilimiting);
            erase(xF, ilimiting);
            AF = cols(A, F);
            AL = cols(A, L);
        }
    };

    // The function that computes the current error norms
    auto update_errors = [&]()
    {
        // Calculate the optimality, feasibility and centrality errors
//        errorf = norminf(gF - tr(AF)*y);
//        errorf = norminf(gF/(tr(AF)*y) - ones(gF.size()));
        errorf = 0.0;
        for(unsigned i = 0; i < F.size(); ++i)
            if(xF[i] > l[F[i]] && errorf < std::abs(gF[i] - dot(AF.col(i), y)))
                errorf = std::abs(gF[i] - dot(AF.col(i), y));
        errorh = norminf(h);

        // Calculate the maximum error
        error = std::max(errorf, errorh);
        result.error = error;
    };

    auto converged = [&]()
    {
        if(error < options.tolerance)
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
        compute_newton_step();
        if(compute_newton_step_failed())
            break;
        update_iterates();
        update_state();
        if(update_state_failed())
            break;
        update_errors();
        output_state();
    } while(!converged());

    outputter.outputHeader();

    // Finish timing the calculation
    result.time = elapsed(begin);

    return result;
}

OptimumSolverActNewton::OptimumSolverActNewton()
: pimpl(new Impl())
{}

OptimumSolverActNewton::OptimumSolverActNewton(const OptimumSolverActNewton& other)
: pimpl(new Impl(*other.pimpl))
{}

OptimumSolverActNewton::~OptimumSolverActNewton()
{}

auto OptimumSolverActNewton::operator=(OptimumSolverActNewton other) -> OptimumSolverActNewton&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto OptimumSolverActNewton::solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult
{
    return pimpl->solve(problem, state, {});
}

auto OptimumSolverActNewton::solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult
{
    return pimpl->solve(problem, state, options);
}

auto OptimumSolverActNewton::dxdp(const Vector& dgdp, const Vector& dbdp) -> Vector
{
    RuntimeError("Could not calculate the sensitivity of the optimal solution with respect to parameters.",
        "The method OptimumSolverActNewton::dxdp has not been implemented yet.");
    return {};
}

auto OptimumSolverActNewton::clone() const -> OptimumSolverBase*
{
    return new OptimumSolverActNewton(*this);
}

} // namespace Reaktoro
