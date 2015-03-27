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
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

/// The function signature of the right-hand side function of a system of ordinary differential equations.
using ODEFunction = std::function<int(double, const Vector&, Vector&)>;

/// The function signature of the Jacobian of the right-hand side function of a system of ordinary differential equations.
using ODEJacobian = std::function<int(double, const Vector&, Matrix&)>;

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

    /// The initial step size to be used in the integration.
    /// An estimation is made if its value is zero.
    double initial_step = 0.0;

    /// The lower bound on the magnitude of the step size.
    double min_step = 0.0;

    /// The upper bound on the magnitude of the step size.
    double max_step = 0.0;

    /// The scalar relative error tolerance.
    double reltol = 1.0e-6;

    /// The scalar absolute error tolerance.
    double abstol = 1.0e-12;

    /// The vector of absolute error tolerances for each component.
    Vector abstols;

    /// The maximum number of steps to be taken by the solver in its attempt to reach the next output time.
    unsigned max_num_steps = 5000;
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
    auto function(double t, const Vector& y, Vector& f) const -> int;

    /// Evaluate the Jacobian of the right-hand side function of the system of ordinary differential equations.
    /// @param t The time variable of the function
    /// @param y The y-variables of the function
    auto jacobian(double t, const Vector& y, Matrix& J) const -> int;

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
    ODESolver(const ODESolver& other);

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
    auto initialize(double tstart, const Vector& y) -> void;

    /// Integrate the ODE performing a single step.
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param[in,out] y The current variables as input, the new current variables as output
    auto integrate(double& t, Vector& y) -> void;

    /// Integrate the ODE performing a single step not going over a given time.
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param[in,out] y The current variables as input, the new current variables as output
    /// @param tfinal The final time that the integration must satisfy
    auto integrate(double& t, Vector& y, double tfinal) -> void;

    /// Solve the ODE equations from a given start time to a final one.
    /// @param tstart The start time for the integration
    /// @param tfinal The final time for the integration
    /// @param[in,out] y The initial value as input, the final value as output
    auto solve(double tstart, double tfinal, Vector& y) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
