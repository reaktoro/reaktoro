// Reaktor is a C++ library for computational reaction modelling.
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

#include "AlgorithmIpnewton.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Common/TimeUtils.hpp>
#include <Reaktor/Optimization/AlgorithmUtils.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>
#include <Reaktor/Optimization/OptimumOptions.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {
namespace {

auto checkInfeasibilityError(const OptimumProblem& problem, const OptimumResult& result) -> void
{
    const auto& lower = problem.lowerBounds();
    const auto& upper = problem.upperBounds();
    if(result.solution.x.size() != problem.numVariables())
        error("Cannot proceed with the minimization.", "Uninitialized primal solution `x`");
    if(result.solution.y.size() != problem.numConstraints())
        error("Cannot proceed with the minimization.", "Uninitialized dual solution `y`");
    if(result.solution.zl.size() != problem.numVariables())
        error("Cannot proceed with the minimization.", "Uninitialized dual solution `zl`");
    if(result.solution.zu.size() != problem.numVariables())
        error("Cannot proceed with the minimization.", "Uninitialized dual solution `zu`");
    if(lower.size() and min(result.solution.x - lower) <= 0.0)
        error("Cannot proceed with the minimization.", "At least one variable in `x` is below its lower bound.");
    if(upper.size() and min(upper - result.solution.x) <= 0.0)
        error("Cannot proceed with the minimization.", "At least one variable in `x` is above its upper bound.");
    if(min(result.solution.zl) < 0.0)
        error("Cannot proceed with the minimization.", "At least one component in `zl` is negative.");
    if(min(result.solution.zu) < 0.0)
        error("Cannot proceed with the minimization.", "At least one component in `zu` is negative.");
}

} // namespace

auto ipnewton(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    Time begin = time();

    const auto& n          = problem.numVariables();
    const auto& m          = problem.numConstraints();
    const auto& objective  = problem.objective();
    const auto& constraint = problem.constraint();
    const auto& lower      = problem.lowerBounds();
    const auto& upper      = problem.upperBounds();
    const auto& tolerance  = options.tolerance;
    const auto& mu         = options.ipnewton.mu;
    const auto& mux        = options.ipnewton.mux;
    const auto& tau        = options.ipnewton.tau;

    Vector& x  = result.solution.x;
    Vector& y  = result.solution.y;
    Vector& zl = result.solution.zl;
    Vector& zu = result.solution.zu;

    const bool has_lower_bounds = lower.size();
    const bool has_upper_bounds = upper.size();

    if(has_lower_bounds) x = max(x, lower + mux*mu*ones(n));
    if(has_upper_bounds) x = max(x, upper - mux*mu*ones(n));

    zl = has_lower_bounds ? Vector(mu/(x - lower).array()) : zeros(n);
    zu = has_upper_bounds ? Vector(mu/(upper - x).array()) : zeros(n);

    checkInfeasibilityError(problem, result);

    ObjectiveResult f;
    ConstraintResult h;

    OptimumStatistics statistics;

    Outputter outputter;

    Vector dx, dy, dzl, dzu;

    Vector dl, du, d;

    auto eval_objective = [&]()
    {
        Time begin = time();
        f = objective(x);
        Time end = time();
        statistics.time_objective_evals += elapsed(end, begin);
    };

    auto eval_contraint = [&]()
    {
        Time begin = time();
        h = constraint(x);
        Time end = time();
        statistics.time_constraint_evals += elapsed(end, begin);
    };

    eval_objective();
    eval_contraint();

    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;

    const SaddlePointOptions& saddle_point_options = options.ipnewton.saddle_point;

    if(options.output.active)
    {
        outputter.setOptions(options.output);

        outputter.addEntry("iter");
        outputter.addEntries("x", n);
        outputter.addEntries("y", m);
        outputter.addEntries("z", n);
        outputter.addEntry("f(x)");
        outputter.addEntry("h(x)");
        outputter.addEntry("errorf");
        outputter.addEntry("errorh");
        outputter.addEntry("errorl");
        outputter.addEntry("erroru");
        outputter.addEntry("error");
        outputter.addEntry("alpha");
        outputter.addEntry("alphaxl");
        outputter.addEntry("alphaxu");
        outputter.addEntry("alphazl");
        outputter.addEntry("alphazu");

        outputter.outputHeader();
        outputter.addValue(statistics.num_iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(zl);
        outputter.addValue(f.func);
        outputter.addValue(norminf(h.func));
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.addValue("---");
        outputter.outputState();
    }

    do
    {
        if(options.ipnewton.scaling)
        {
            dl = has_lower_bounds ? sqrt(x - lower).eval() : ones(n);
            du = has_upper_bounds ? sqrt(upper - x).eval() : ones(n);
            d  = dl % du;

            const auto D = diag(d);

            saddle_point_problem.H = D*f.hessian*D;
            saddle_point_problem.A = h.grad*D;
            saddle_point_problem.f = -d % f.grad + D*h.grad.transpose()*y;
            saddle_point_problem.g = -h.func;

            if(has_lower_bounds and has_upper_bounds) {
                diagonal(saddle_point_problem.H) += (upper - x)%zl + (x - lower)%zu;
                saddle_point_problem.f += mu*(du/dl - dl/du); }
            if(has_lower_bounds and not has_upper_bounds) {
                diagonal(saddle_point_problem.H) += zl;
                saddle_point_problem.f += mu/dl; }
            if(has_upper_bounds and not has_lower_bounds) {
                diagonal(saddle_point_problem.H) += zu;
                saddle_point_problem.f -= mu/du; }

            solve(saddle_point_problem, saddle_point_result, saddle_point_options);

            statistics.time_linear_system_solutions += saddle_point_result.statistics.time;

            dx = d % saddle_point_result.solution.x;
            dy = saddle_point_result.solution.y;
            if(has_lower_bounds) dzl = (mu - zl%dx)/(x - lower) - zl;
            if(has_upper_bounds) dzu = (mu + zu%dx)/(upper - x) - zu;
        }
        else
        {
            saddle_point_problem.H = f.hessian;
            saddle_point_problem.H.diagonal() += zl/x;
            saddle_point_problem.A = h.grad;
            saddle_point_problem.f = -(f.grad - h.grad.transpose()*y - mu/x);
            saddle_point_problem.g = -h.func;

            solve(saddle_point_problem, saddle_point_result, saddle_point_options);

            dx = saddle_point_result.solution.x;
            dy = saddle_point_result.solution.y;
            dzl = (mu - zl%dx)/x - zl;
        }

        const double alphaxl = has_lower_bounds ? fractionToTheBoundary(x - lower, dx, tau) : 1.0;
        const double alphaxu = has_upper_bounds ? fractionToTheBoundary(upper - x, dx, tau) : 1.0;
        const double alphazl = has_lower_bounds ? fractionToTheBoundary(zl, dzl, tau) : 1.0;
        const double alphazu = has_upper_bounds ? fractionToTheBoundary(zl, dzu, tau) : 1.0;

        double alpha;

        if(options.ipnewton.uniform_newton_step)
        {
            alpha = std::min({alphaxl, alphaxu, alphazl, alphazu});

            x += alpha * dx;
            y += alpha * dy;
            if(has_lower_bounds) zl += alpha * dzl;
            if(has_upper_bounds) zu += alpha * dzu;
        }
        else
        {
            alpha = std::min({alphaxl, alphaxu});

            x += alpha * dx;
            y += dy;
            if(has_lower_bounds) zl += alphazl * dzl;
            if(has_upper_bounds) zu += alphazu * dzu;
        }

        eval_objective();
        eval_contraint();
//        f = objective(x);
//        h = constraint(x);

        // Calculate the optimality, feasibility and centrality errors
        const double errorf = norminf(f.grad - h.grad.transpose()*y - zl + zu);
        const double errorh = norminf(h.func);
        const double errorl = has_lower_bounds ? norminf((x - lower) % zl - mu) : 0.0;
        const double erroru = has_upper_bounds ? norminf((upper - x) % zu - mu) : 0.0;

        // Calculate the maximum error
        statistics.error = std::max({errorf, errorh, errorl, erroru});

        ++statistics.num_iterations;

        if(options.output.active)
        {
            outputter.addValue(statistics.num_iterations);
            outputter.addValues(x);
            outputter.addValues(y);
            outputter.addValues(zl);
            outputter.addValue(f.func);
            outputter.addValue(norminf(h.func));
            outputter.addValue(errorf);
            outputter.addValue(errorh);
            outputter.addValue(errorl);
            outputter.addValue(erroru);
            outputter.addValue(statistics.error);
            outputter.addValue(alpha);
            outputter.addValue(alphaxl);
            outputter.addValue(alphaxu);
            outputter.addValue(alphazl);
            outputter.addValue(alphazu);
            outputter.outputState();
        }

    } while(statistics.error > tolerance and statistics.num_iterations < options.max_iterations);

    outputter.outputHeader();

    if(statistics.num_iterations < options.max_iterations)
        statistics.converged = true;

    result.statistics = statistics;

    Time end = time();

    result.statistics.time = elapsed(end, begin);
}
