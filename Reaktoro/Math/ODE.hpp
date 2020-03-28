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

#pragma once

// C++ includes
#include <functional>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The function signature of the right-hand side function of a system of ordinary differential equations.
using ODEFunction = std::function<int(double, VectorXdConstRef, VectorXdRef)>;

/// The function signature of the Jacobian of the right-hand side function of a system of ordinary differential equations.
using ODEJacobian = std::function<int(double, VectorXdConstRef, MatrixXdRef)>;

/// The linear multistep method to be used in ODESolver.
enum class ODEStepMode { Adams, BDF };

/// The type of nonlinear solver iteration to be used in ODESolver.
enum class ODEIterationMode { Functional, Newton };

/// A struct that defines the options for the ODESolver.
/// @see ODESolver, ODEProblem
struct ODEOptions
{
    /// The linear multistep method used in the integration.
    ODEStepMode step = ODEStepMode::BDF;

    /// The type of nonlinear solver iteration used in the integration.
    ODEIterationMode iteration = ODEIterationMode::Newton;

    /// The flag that enables the STAbility Limit Detection (STALD) algorithm.
    /// The STALD algorithm should be used when BDF method does not progress well,
    /// which can happen when the current BDF order is above 2. Using the STALD
    /// algorithm helps fixing this by reducing the current BDF order, and thus
    /// increasing its stability.
    bool stability_limit_detection = false;

    /// The initial step size to be used in the integration.
    /// An estimation is made if its value is zero.
    double initial_step = 0.0;

    /// The value of the independent variable t past which the solution is not to proceed.
    /// No stop time is used if its value is zero.
    double stop_time = 0.0;

    /// The lower bound on the magnitude of the step size.
    double min_step = 0.0;

    /// The upper bound on the magnitude of the step size.
    double max_step = 0.0;

    /// The scalar relative error tolerance.
    double reltol = 1e-4;

    /// The scalar absolute error tolerance.
    double abstol = 1.0e-6;

    /// The maximum order for the BDF integration scheme.
    /// The order of the BDF method can be any integer from 1 to 5 inclusive.
    unsigned max_order_bdf = 5;

    /// The maximum order for the Adams integration scheme.
    /// The order of the Adams method can be any integer from 1 to 12 inclusive.
    unsigned max_order_adams = 5;

    /// The maximum allowed number of steps before reaching the final time.
    unsigned max_num_steps = 500;

    /// The maximum number of warnings for `t + h = t`, with `h` being too small compared to `t`.
    unsigned max_hnil_warnings = 10;

    /// The maximum number of error test failures.
    unsigned max_num_error_test_failures = 20;

    /// The maximum number of nonlinear iterations.
    unsigned max_num_nonlinear_iterations = 3;

    /// The maximum number of convergence failures.
    unsigned max_num_convergence_failures = 10;

    /// The coefficient in the nonlinear convergence test.
    double nonlinear_convergence_coefficient = 0.1;

    /// The vector of absolute error tolerances for each component.
    VectorXd abstols;
};

/// A class that defines a system of ordinary differential equations (ODE) problem.
/// @see ODESolver, ODEOptions
class ODEProblem
{
public:
    /// Construct a default ODEProblem instance
    ODEProblem();

    /// Construct a copy of an ODEProblem instance.
    ODEProblem(const ODEProblem& other);

    /// Destroy this ODEProblem instance.
    virtual ~ODEProblem();

    /// Assign a ODEProblem instance to this instance
    auto operator=(ODEProblem other) -> ODEProblem&;

    /// Set the number of ordinary differential equations
    auto setNumEquations(unsigned num) -> void;

    /// Set the right-hand side function of the system of ordinary differential equations
    auto setFunction(const ODEFunction& f) -> void;

    /// Set the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto setJacobian(const ODEJacobian& J) -> void;

    /// Return true if the problem has bee initialized.
    auto initialized() const -> bool;

    /// Return the number of ordinary differential equations
    auto numEquations() const -> unsigned;

    /// Return the right-hand side function of the system of ordinary differential equations
    auto function() const -> const ODEFunction&;

    /// Return the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto jacobian() const -> const ODEJacobian&;

    /// Evaluate the right-hand side function of the system of ordinary differential equations.
    /// @param t The time variable of the function
    /// @param y The y-variables of the function
    /// @param[out] f The result of the function evaluation.
    /// @return Return 0 if successful, any other number otherwise.
    auto function(double t, VectorXdConstRef y, VectorXdRef f) const -> int;

    /// Evaluate the Jacobian of the right-hand side function of the system of ordinary differential equations.
    /// @param t The time variable of the function
    /// @param y The y-variables of the function
    /// @param[out] J The result of the Jacobian evaluation.
    /// @return Return 0 if successful, any other number otherwise.
    auto jacobian(double t, VectorXdConstRef y, MatrixXdRef J) const -> int;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// A wrapper class for `CVODE`, a library for solving ordinary differential equations.
/// @see ODEProblem, ODEOptions
class ODESolver
{
public:
    /// Construct a default ODESolver instance.
    ODESolver();

    /// Construct a copy of an ODESolver instance.
    ODESolver(const ODESolver& other) = delete;

    /// Destroy this ODESolver instance.
    virtual ~ODESolver();

    /// Assign a ODESolver instance to this instance
    auto operator=(ODESolver other) -> ODESolver&;

    /// Set the options for the ODE solver.
    /// @see ODEOptions
    auto setOptions(const ODEOptions& options) -> void;

    /// Set the ODE problem.
    /// @see ODEProblem
    auto setProblem(const ODEProblem& problem) -> void;

    /// Initializes the ODE solver.
    /// This method should be invoked whenever the user intends to make a call to `ODESolver::integrate`.
    /// @param tstart The start time of the integration.
    /// @param y The initial values of the variables
    auto initialize(double tstart, VectorXdConstRef y) -> void;

    /// Integrate the ODE performing a single step.
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param[in,out] y The current variables as input, the new current variables as output
    auto integrate(double& t, VectorXdRef y) -> void;

    /// Integrate the ODE performing a single step not going over a given time.
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param[in,out] y The current variables as input, the new current variables as output
    /// @param tfinal The final time that the integration must satisfy
    auto integrate(double& t, VectorXdRef y, double tfinal) -> void;

    /// Solve the ODE equations from a given start time to a final one.
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param dt The value of the time step
    /// @param[in,out] y The current variables as input, the new current variables as output
    auto solve(double& t, double dt, VectorXdRef y) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
