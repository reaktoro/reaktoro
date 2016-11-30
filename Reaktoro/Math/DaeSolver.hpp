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

#pragma once

// C++ includes
#include <functional>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The function signature of the right-hand side function of a system of ordinary differential equations.
using DaeFunction = std::function<int(double t, const Vector& u, const Vector& udot, Vector& F)>;

/// The function signature of the Jacobian of the right-hand side function of a system of ordinary differential equations.
using DaeJacobian = std::function<int(double t, const Vector& u, const Vector& udot, Matrix& J)>;

/// A struct that defines the options for the DaeSolver.
/// @see DaeSolver, DaeProblem
struct DaeOptions
{
};

/// A class that defines a system of differential-algebraic equations (DAE) problem.
/// @see DaeSolver, DaeOptions
class DaeProblem
{
public:
    /// Construct a default DaeProblem instance
    DaeProblem();

    /// Construct a copy of an DaeProblem instance.
    DaeProblem(const DaeProblem& other);

    /// Destroy this DaeProblem instance.
    virtual ~DaeProblem();

    /// Assign a DaeProblem instance to this instance
    auto operator=(DaeProblem other) -> DaeProblem&;

    /// Set the number of ordinary differential equations
    auto setNumEquations(unsigned num) -> void;

    /// Set the right-hand side function of the system of ordinary differential equations
    auto setFunction(const DaeFunction& f) -> void;

    /// Set the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto setJacobian(const DaeJacobian& J) -> void;

    /// Return true if the problem has bee initialized.
    auto initialized() const -> bool;

    /// Return the number of ordinary differential equations
    auto numEquations() const -> unsigned;

    /// Return the right-hand side function of the system of ordinary differential equations
    auto function() const -> const DaeFunction&;

    /// Return the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto jacobian() const -> const DaeJacobian&;

    /// Evaluate the right-hand side function of the system of ordinary differential equations.
    /// @param t The independent progress variable.
    /// @param y The values of the y variables.
    /// @param ydot The values of the first-order derivatives of the y variables.
    /// @param[out] f The result of the function evaluation.
    /// @return Return 0 if successful, any other number otherwise.
    auto function(double t, const Vector& y, const Vector& ydot, Vector& f) const -> int;

    /// Evaluate the Jacobian of the right-hand side function of the system of ordinary differential equations.
    /// @param t The independent progress variable.
    /// @param y The values of the y variables.
    /// @param ydot The values of the first-order derivatives of the y variables.
    /// @param[out] J The result of the Jacobian evaluation.
    /// @return Return 0 if successful, any other number otherwise.
    auto jacobian(double t, const Vector& y, const Vector& ydot, Matrix& J) const -> int;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// A class for solving differential-algebraic equations (DAE).
/// @see DaeProblem, DaeOptions
class DaeSolver
{
public:
    /// Construct a default DaeSolver instance.
    DaeSolver();

    /// Construct a copy of an DaeSolver instance.
    DaeSolver(const DaeSolver& other);

    /// Destroy this DaeSolver instance.
    virtual ~DaeSolver();

    /// Assign a DaeSolver instance to this instance
    auto operator=(DaeSolver other) -> DaeSolver&;

    /// Set the options for the ODE solver.
    /// @see DaeOptions
    auto setOptions(const DaeOptions& options) -> void;

    /// Set the ODE problem.
    /// @see DaeProblem
    auto setProblem(const DaeProblem& problem) -> void;

    /// Initializes the ODE solver.
    /// This method should be invoked whenever the user intends to make a call to `DaeSolver::integrate`.
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
    /// @param[in,out] t The current time of the integration as input, the new current time as output
    /// @param dt The value of the time step
    /// @param[in,out] y The current variables as input, the new current variables as output
    auto solve(double& t, double dt, Vector& y) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
