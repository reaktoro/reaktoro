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

namespace Reaktor {

/// The function signature of the right-hand side function of
/// a system of ordinary differential equations.
typedef std::function<int(double, const Vector&, Vector&)> ODEFunction;

/// The function signature of the Jacobian of the right-hand side function
/// of a system of ordinary differential equations.
typedef std::function<int(double, const Vector&, Matrix&)> ODEJacobian;

/// A class that defines a system of ordinary differential equations (ODE) problem.
class ODEProblem
{
public:
    /// Convenient type definition for ODEFunction
    typedef ODEFunction Function;

    /// Convenient type definition for @ref ODEJacobian
    typedef ODEJacobian Jacobian;

    /// Construct a default @ref ODEProblem instance
    ODEProblem();

    /// Set the number of ordinary differential equations
    auto setNumEquations(unsigned num) -> void;

    /// Set the right-hand side function of the system of ordinary differential equations
    auto setFunction(const Function& f) -> void;

    /// Set the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto setJacobian(const Jacobian& J) -> void;

    /// Return the number of ordinary differential equations
    auto numEquations() const -> unsigned;

    /// Return the right-hand side function of the system of ordinary differential equations
    auto function() const -> const Function&;

    /// Return the Jacobian of the right-hand side function of the system of ordinary differential equations
    auto jacobian() const -> const Jacobian&;

    /// Evaluate the right-hand side function of the system of ordinary differential equations.
    /// @param t The time variable of the function
    /// @param y The y-variables of the function
    auto function(double t, const Vector& y, Vector& f) const -> int;

    /// Evaluate the Jacobian of the right-hand side function of the system of ordinary differential equations.
    /// @param t The time variable of the function
    /// @param y The y-variables of the function
    auto jacobian(double t, const Vector& y, Matrix& J) const -> int;

private:
    /// The number of ordinary differential equations
    unsigned num_equations = 0;

    /// The right-hand side function of the system of ordinary differential equations
    Function _function;

    /// The Jacobian of the right-hand side function of the system of ordinary differential equations
    Jacobian _jacobian;
};

} // namespace Reaktor
