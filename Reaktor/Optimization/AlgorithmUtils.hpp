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

#pragma once

// C++ includes
#include <functional>
#include <list>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Cube.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/// A type that describes the result of the evaluation of a objective function
/// @see ObjectiveFunction
struct ObjectiveResult
{
    /// The result of the evaluation of the objective function
    double func;

    /// The gradient result of the evaluation of the objective function
    Vector grad;

    /// The Hessian result of the evaluation of the objective function
    Matrix hessian;
};

/// A type that describes the result of the evaluation of a constraint function
/// @see ConstraintFunction
struct ConstraintResult
{
    /// The result of the evaluation of the constraint function
    Vector func;

    /// The gradient result of the evaluation of the constraint function
    Matrix grad;

    /// The Hessian result of the evaluation of the constraint function
    Cube hessian;
};

/// A type that describes the functional signature of an objective function
/// @see ObjectiveResult
typedef std::function<
    ObjectiveResult(const Vector&)>
        ObjectiveFunction;

/// A type that describes the functional signature of a constraint function
/// @see ConstraintResult
typedef std::function<
    ConstraintResult(const Vector&)>
        ConstraintFunction;

/// A type that describes the definition of an optimization problem
class OptimumProblem
{
public:
    /// Construct a OptimumProblem instance
    /// @param n The number of primal variables
    /// @param m The number of equality constraints
    OptimumProblem(unsigned n, unsigned m);

    /// Set the objective function of the optimization problem
    /// @param f The objective function of the optimization problem
    auto setObjective(const ObjectiveFunction& objective) -> void;

    /// Set the equality constraint function of the optimization problem
    /// @param A The coefficient matrix of the equality constraints
    /// @param b The right-hand side vector of the equality constraints
    auto setConstraint(const ConstraintFunction& constraint) -> void;

    /// Set the lower bounds of the optimization problem
    /// @param l The lower bounds of the primal variables
    auto setLowerBounds(const Vector& l) -> void;

    /// Set the lower bounds of the optimization problem
    /// @param u The upper bounds of the primal variables
    auto setUpperBounds(const Vector& u) -> void;

    /// Get the number of variables in the optimization problem
    auto numVariables() const -> unsigned;

    /// Get the number of equality constraints in the optimization problem
    auto numConstraints() const -> unsigned;

    /// Get the objective function of the optimization problem
    auto objective() const -> const ObjectiveFunction&;

    /// Get the equality constraint function of the optimization problem
    auto constraint() const -> const ConstraintFunction&;

    /// Get the lower bounds of the optimization problem
    auto lowerBounds() const -> const Vector&;

    /// Get the upper bounds of the optimization problem
    auto upperBounds() const -> const Vector&;

private:
    /// The number of primal variables
    unsigned n;

    /// The number of equality constraints
    unsigned m;

    /// The objective function to be minimized in the optimization problem
    ObjectiveFunction f;

    /// The coefficient matrix of the equality constraints in the optimization problem
    ConstraintFunction h;

    /// The lower bound vector of the inequality constraints in the optimization problem
    Vector l;

    /// The upper bound vector of the inequality constraints in the optimization problem
    Vector u;
};

/// A type that describes the solution of an optimization problem
struct OptimumSolution
{
    /// The primal solution of the optimization problem
    Vector x;

    /// The dual solution of the optimization problem with respect to the equality constraints
    Vector y;

    /// The dual solution of the optimization problem with respect to the lower bound constraints
    Vector zl;

    /// The dual solution of the optimization problem with respect to the upper bound constraints
    Vector zu;
};

/// A type that describes the result of an optimization calculation
struct OptimumStatistics
{
    /// The flag that indicates if the optimization calculation converged
    bool converged;

    /// The number of iterations in the optimization calculation
    unsigned num_iterations;

    /// The number of evaluations of the objective function in the optimization calculation
    unsigned num_objective_evals;

    /// The convergence rate of the optimization calculation near the solution
    double convergence_rate;

    /// The final residual error of the optimization calculation
    double error;
};

/// A type that contains important internal state for the optimization calculation.
/// @warning The internal state of this class must not be changed outside the optimization algorithms.
struct OptimumInternal
{

};

/// A type that describes the result of an optimization calculation
struct OptimumResult
{
    /// The solution of the optimization calculation
    OptimumSolution solution;

    /// The statistics of the optimization calculation
    OptimumStatistics statistics;

    /// The internal state of the optimization calculation
    OptimumInternal internal;
};

struct IpoptParams
{
    double mu              = 1.0e-8;
    double delta           = 1.0;
    double eta_phi         = 1.0e-4;
    double gamma_alpha     = 0.05;
    double gamma_phi       = 1.0e-5;
    double gamma_theta     = 1.0e-5;
    double kappa_epsilon   = 10.0;
    double kappa_mu        = 0.2;
    double kappa_sigma     = 100;
    double kappa_soc       = 0.99;
    double ksi_phi         = 1.0e-15;
    double s_phi           = 2.3;
    double s_theta         = 1.1;
    double tau_min         = 0.99;
    double theta_mu        = 1.5;
    unsigned max_iters_soc = 4;
};

/// A type that describes
struct OptimumOptions
{
    /// The residual tolerance in the optimization calculation
    double tolerance = 1.0e-8;

    /// The interior-point method perturbation in the optimization calculation
    double mu = 1.0e-8;

    IpoptParams ipopt;
};

/// A type that describes the entry of an optimization filter
typedef std::vector<double> FilterEntry;

/// A type that describes an optimization filter
typedef std::list<FilterEntry> Filter;

/// Check if a filter entry is dominated by another.
/// Check if entry `a` is dominated by `b`, that is, if `a` > `b` componentwise.
/// @param a The filter entry `a`
/// @param b The filter entry `b`
auto dominated(const FilterEntry& a, const FilterEntry& b) -> bool;

/// Check if an entry is acceptable to a filter
/// @param entry The entry to be checked
/// @param filter The filter where the entry is checked
auto acceptable(const Filter& filter, const FilterEntry& entry) -> bool;

/// Extend a filter with a new entry.
/// This method not only insert a new entry to the filter,
/// but also removes all existent dominated entries with
/// the insertion of the new one.
/// @param filter The filter where the entry is inserted
/// @param entry The entry to be inserted in the filter
auto extend(Filter& filter, const FilterEntry& entry) -> void;

/// Compute the largest step length @f$\alpha@f$ such that
/// @f$\mathbf{p} + \alpha\Delta\mathbf{p}@f$ is on the
/// lower bound @f$\mathbf{x}_l=\mathbf{0}@f$.
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
auto largestStep(const Vector& p, const Vector& dp) -> double;

/// Compute the fraction-to-the-boundary step length given by:
/// @f[\alpha_{\mathrm{max}}=\max\{\alpha\in(0,1]:\mathbf{p}+\alpha\Delta\mathbf{p}\geq(1-\tau)\mathbf{p}\}@f.]
/// @param p The point @f$\mathbf{p}@f$
/// @param dp The step @f$\Delta\mathbf{p}@f$
/// @param tau The fraction-to-the-boundary parameter @f$\tau@f$
auto fractionToTheBoundary(const Vector& p, const Vector& dp, double tau) -> double;

/// Check if a float number is less than another by a base value.
/// This method is particularly useful when comparing two float numbers
/// where round-off errors can compromise the checking.
/// The following is used for the comparison:
/// @f[a < b + 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto lessThan(double a, double b, double baseval) -> bool;

/// Check if a float number is greater than another by a base value.
/// This method is particularly useful when comparing two float numbers
/// where round-off errors can compromise the checking.
/// The following is used for the comparison:
/// @f[a > b - 10\epsilon \mathrm{baseval}@f],
/// where @f$\epsilon@f$ is the machine double precision.
/// @param a The left-hand side value
/// @param b The right-hand side value
/// @param baseval The base value for the comparison
auto greaterThan(double a, double b, double baseval) -> bool;

} // namespace Reaktor
