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

// Armadillo includes
#include <armadillo>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/Outputter.hpp>
#include <Reaktor/Optimization/AlgorithmUtils.hpp>
#include <Reaktor/Optimization/SaddlePointUtils.hpp>

namespace Reaktor {

auto ipnewton(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
{
    const auto& n          = problem.numVariables();
    const auto& m          = problem.numConstraints();
    const auto& objective  = problem.objective();
    const auto& constraint = problem.constraint();
    const auto& lower      = problem.lowerBounds();
    const auto& upper      = problem.upperBounds();
    const auto& tolerance  = options.tolerance;
    const auto& mu         = options.ipopt.mu;
    const auto& tau        = options.ipopt.tau_min;

    Vector& x  = result.solution.x;
    Vector& y  = result.solution.y;
    Vector& zl = result.solution.zl;
    Vector& zu = result.solution.zu;

    const bool has_lower_bounds = not lower.empty();
    const bool has_upper_bounds = not upper.empty();

    if(not has_lower_bounds) zl = arma::zeros(n);
    if(not has_upper_bounds) zu = arma::zeros(n);

    ObjectiveResult f;
    ConstraintResult h;

    Vector dx, dy, dzl, dzu;

    f = objective(x);
    h = constraint(x);

    SaddlePointProblem saddle_point_problem;
    SaddlePointResult saddle_point_result;

    OptimumStatistics statistics;

    Outputter outputter;

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
    outputter.addValue(arma::norm(h.func, "inf"));
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.addValue("---");
    outputter.outputState();

    do
    {
        const Vector dl = has_lower_bounds ? arma::sqrt(x - lower).eval() : arma::ones(n);
        const Vector du = has_upper_bounds ? arma::sqrt(upper - x).eval() : arma::ones(n);
        const Vector d  = dl % du;

        const auto D = arma::diagmat(d);

        saddle_point_problem.H = D*f.hessian*D;
        saddle_point_problem.A = h.grad*D;
        saddle_point_problem.f = -D*f.grad + D*h.grad.t()*y;
        saddle_point_problem.g = -h.func;

        if(has_lower_bounds and has_upper_bounds) {
            saddle_point_problem.H.diag() += (upper - x)%zl + (x - lower)%zu;
            saddle_point_problem.f += mu*(du/dl - dl/du); }
        if(has_lower_bounds and not has_upper_bounds) {
            saddle_point_problem.H.diag() += zl;
            saddle_point_problem.f += mu/dl; }
        if(has_upper_bounds and not has_lower_bounds) {
            saddle_point_problem.H.diag() += zu;
            saddle_point_problem.f -= mu/du; }

        solveNullspace(saddle_point_problem, saddle_point_result);

        dx = d % saddle_point_result.solution.x;
        dy = saddle_point_result.solution.y;
        if(has_lower_bounds) dzl = (mu - zl%dx)/(x - lower) - zl;
        if(has_upper_bounds) dzu = (mu + zu%dx)/(upper - x) - zu;

        const double alphaxl = has_lower_bounds ? fractionToTheBoundary(x - lower, dx, tau) : 1.0;
        const double alphaxu = has_upper_bounds ? fractionToTheBoundary(upper - x, dx, tau) : 1.0;
        const double alphazl = has_lower_bounds ? fractionToTheBoundary(zl, dzl, tau) : 1.0;
        const double alphazu = has_upper_bounds ? fractionToTheBoundary(zl, dzu, tau) : 1.0;

        const double alpha = std::min({alphaxl, alphaxu, alphazl, alphazu});

        x += alpha * dx;
        y += alpha * dy;
        if(has_lower_bounds) zl += alpha * dzl;
        if(has_upper_bounds) zu += alpha * dzu;

        f = objective(x);
        h = constraint(x);

        // Calculate the optimality, feasibility and centrality errors
        const double errorf = arma::norm(f.grad - h.grad.t()*y - zl + zu, "inf");
        const double errorh = arma::norm(h.func, "inf");
        const double errorl = has_lower_bounds ? arma::norm((x - lower) % zl - mu, "inf") : 0.0;
        const double erroru = has_upper_bounds ? arma::norm((upper - x) % zu - mu, "inf") : 0.0;

        // Calculate the maximum error
        statistics.error = std::max({errorf, errorh, errorl, erroru});

        ++statistics.num_iterations;

        outputter.addValue(statistics.num_iterations);
        outputter.addValues(x);
        outputter.addValues(y);
        outputter.addValues(zl);
        outputter.addValue(f.func);
        outputter.addValue(arma::norm(h.func, "inf"));
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

    } while(statistics.error > tolerance and statistics.num_iterations < options.max_iterations);

    outputter.outputHeader();

    if(statistics.num_iterations < options.max_iterations)
        statistics.converged = true;

    result.statistics = statistics;
}

//auto ipnewton(const OptimumProblem& problem, OptimumResult& result, const OptimumOptions& options) -> void
//{
//    const auto& objective  = problem.objective();
//    const auto& constraint = problem.constraint();
//    const auto& tolerance  = options.tolerance;
//    const auto& mu         = options.ipopt.mu;
//    const auto& tau        = options.ipopt.tau_min;
//
//    Vector& x  = result.solution.x;
//    Vector& y  = result.solution.y;
//    Vector& zl = result.solution.zl;
//
//    const Vector& lower = problem.lowerBounds();
//
//    ObjectiveResult f, phi;
//    ConstraintResult h;
//
//    Vector dx, dy, dz;
//
//    f = objective(x);
//    h = constraint(x);
//
//    SaddlePointProblem saddle_point_problem;
//    SaddlePointResult saddle_point_result;
//
//    OptimumStatistics statistics;
//
//    Outputter outputter;
//
//    outputter.setOptions(options.output);
//
//    const unsigned n = problem.numVariables();
//    const unsigned m = problem.numConstraints();
//
//    outputter.addEntry("iter");
//    outputter.addEntries("x", n);
//    outputter.addEntries("y", m);
//    outputter.addEntries("z", n);
//    outputter.addEntry("f(x)");
//    outputter.addEntry("h(x)");
//    outputter.addEntry("mu(w)");
//    outputter.addEntry("errorf");
//    outputter.addEntry("errorh");
//    outputter.addEntry("errorc");
//    outputter.addEntry("error");
//    outputter.addEntry("alpha");
//    outputter.addEntry("alphax");
//    outputter.addEntry("alphaz");
//
//    outputter.outputHeader();
//    outputter.addValue(statistics.num_iterations);
//    outputter.addValues(x);
//    outputter.addValues(y);
//    outputter.addValues(zl);
//    outputter.addValue(f.func);
//    outputter.addValue(arma::norm(h.func, "inf"));
//    outputter.addValue(arma::dot((x - lower), zl)/n);
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.addValue("---");
//    outputter.outputState();
//
//    do
//    {
//        phi.func = f.func - mu * arma::sum(arma::log(x - lower));
//        phi.grad = f.grad - mu/(x - lower);
//
//        saddle_point_problem.A = h.grad;
//        saddle_point_problem.H = f.hessian + arma::diagmat(zl/(x - lower));
//        saddle_point_problem.f = -(phi.grad - h.grad.t()*y);
//        saddle_point_problem.g = -h.func;
//
//        solveNullspace(saddle_point_problem, saddle_point_result);
//
//        dx = saddle_point_result.solution.x;
//        dy = saddle_point_result.solution.y;
//        dz = -(zl % dx + (x - lower) % zl - mu)/(x - lower);
//
//        double alphax = fractionToTheBoundary(x - lower, dx, tau);
//        double alphaz = fractionToTheBoundary(zl, dz, tau);
//        double alpha = std::min(alphax, alphaz);
//
//        x  += alpha * dx;
//        y  += alpha * dy;
//        zl += alpha * dz;
//
//        f = objective(x);
//        h = constraint(x);
//
//        // Calculate the optimality, feasibility and centrality errors
//        const double errorf = arma::norm(f.grad - h.grad.t()*y - zl, "inf");
//        const double errorh = arma::norm(h.func, "inf");
//        const double errorc = arma::norm((x - lower) % zl - mu, "inf");
//
//        // Calculate the maximum error
//        statistics.error = std::max({errorf, errorh, errorc});
//
//        ++statistics.num_iterations;
//
//        outputter.addValue(statistics.num_iterations);
//        outputter.addValues(x);
//        outputter.addValues(y);
//        outputter.addValues(zl);
//        outputter.addValue(f.func);
//        outputter.addValue(arma::norm(h.func, "inf"));
//        outputter.addValue(arma::dot((x - lower), zl)/n);
//        outputter.addValue(errorf);
//        outputter.addValue(errorh);
//        outputter.addValue(errorc);
//        outputter.addValue(statistics.error);
//        outputter.addValue(alpha);
//        outputter.addValue(alphax);
//        outputter.addValue(alphaz);
//        outputter.outputState();
//
//    } while(statistics.error > tolerance and statistics.num_iterations < options.max_iterations);
//
//    outputter.outputHeader();
//
//    if(statistics.num_iterations < options.max_iterations)
//        statistics.converged = true;
//
//    result.statistics = statistics;
//}

} // namespace Reaktor
