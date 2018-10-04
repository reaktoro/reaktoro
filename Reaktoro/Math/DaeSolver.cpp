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

#include "DaeSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct DaeProblem::Impl
{
    /// The number of ordinary differential equations
    unsigned num_equations = 0;

    /// The right-hand side function of the system of ordinary differential equations
    DaeFunction ode_function;

    /// The Jacobian of the right-hand side function of the system of ordinary differential equations
    DaeJacobian ode_jacobian;
};

struct DaeSolver::Impl
{
    /// The DAE problem
    DaeProblem problem;

    /// The options for the DAE integration
    DaeOptions options;

    /// The auxiliary vector y
    Vector y;

    /// The auxiliary vector f for the function evaluation
    Vector f;

    /// The auxiliary matrix J for the Jacobian evaluation
    Matrix J;

    /// Construct a default DaeSolver::Impl instance
    Impl()
    {}

    ~Impl()
    {
    }

    /// Initializes the DAE solver
    auto initialize(double tstart, VectorConstRef y) -> void
    {
        // Check if the ordinary differential problem has been initialized
        Assert(problem.initialized(),
            "Cannot proceed with DaeSolver::initialize to initialize the solver.",
            "The provided DaeProblem instance was not properly initialized.");

        // Check if the dimension of 'y' matches the number of equations
        Assert(y.size() == problem.numEquations(),
            "Cannot proceed with DaeSolver::initialize to initialize the solver.",
            "The dimension of the vector parameter `y` does not match the number of equations.");
    }

    /// Integrate the DAE performing a single step.
    auto integrate(double& t, VectorRef y) -> void
    {
    }

    /// Integrate the DAE performing a single step not going over a given time.
    auto integrate(double& t, VectorRef y, double tfinal) -> void
    {
    }

    /// Solve the DAE equations from a given start time to a final one.
    auto solve(double& t, double dt, VectorRef y) -> void
    {
    }
};

DaeProblem::DaeProblem()
: pimpl(new Impl())
{}

DaeProblem::DaeProblem(const DaeProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

DaeProblem::~DaeProblem()
{}

auto DaeProblem::operator=(DaeProblem other) -> DaeProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto DaeProblem::setNumEquations(unsigned num) -> void
{
    pimpl->num_equations = num;
}

auto DaeProblem::setFunction(const DaeFunction& f) -> void
{
    pimpl->ode_function = f;
}

auto DaeProblem::setJacobian(const DaeJacobian& J) -> void
{
    pimpl->ode_jacobian = J;
}

auto DaeProblem::initialized() const -> bool
{
    return numEquations() && function();
}

auto DaeProblem::numEquations() const -> unsigned
{
    return pimpl->num_equations;
}

auto DaeProblem::function() const -> const DaeFunction&
{
    return pimpl->ode_function;
}

auto DaeProblem::jacobian() const -> const DaeJacobian&
{
    return pimpl->ode_jacobian;
}

auto DaeProblem::function(double t, VectorConstRef y, VectorConstRef ydot, VectorRef f) const -> int
{
    return function()(t, y, ydot, f);
}

auto DaeProblem::jacobian(double t, VectorConstRef y, VectorConstRef ydot, MatrixRef J) const -> int
{
    return jacobian()(t, y, ydot, J);
}

DaeSolver::DaeSolver()
: pimpl(new Impl())
{}

DaeSolver::DaeSolver(const DaeSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

DaeSolver::~DaeSolver()
{}

auto DaeSolver::operator=(DaeSolver other) -> DaeSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto DaeSolver::setOptions(const DaeOptions& options) -> void
{
    pimpl->options = options;
}

auto DaeSolver::setProblem(const DaeProblem& problem) -> void
{
    pimpl->problem = problem;
}

auto DaeSolver::initialize(double tstart, VectorConstRef y) -> void
{
    pimpl->initialize(tstart, y);
}

auto DaeSolver::integrate(double& t, VectorRef y) -> void
{
    pimpl->integrate(t, y);
}

auto DaeSolver::integrate(double& t, VectorRef y, double tfinal) -> void
{
    pimpl->integrate(t, y, tfinal);
}

auto DaeSolver::solve(double& t, double dt, VectorRef y) -> void
{
    pimpl->solve(t, dt, y);
}

} // namespace Reaktoro
