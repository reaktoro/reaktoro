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

#include "EquilibriumProblem.hpp"

// Reaktor includes
#include <Reaktor/Common/MatrixUtils.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>

namespace Reaktor {

struct EquilibriumProblem::Impl
{
    /// The reference to the ChemicalSystem instance
    const ChemicalSystem& system;

    /// The reference to the Partition instance
    const Partition& partition;

    /// The temperature for the equilibrium problem (in units of K)
    double T;

    /// The pressure for the equilibrium problem (in units of Pa)
    double P;

    /// The electrical charge for the equilibrium problem (in units of mol)
    double charge;

    /// The amounts of the elements for the equilibrium problem (in units of mol)
    Vector b;

    /// The balance matrix of the system with linearly independent rows
    Matrix A;

    /// The indices of the linearly independent components
    Indices independent_components;

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const ChemicalSystem& system)
    : Impl(system, Partition::allEquilibrium(system))
    {}

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const ChemicalSystem& system, const Partition& partition)
    : system(system), partition(partition),
      T(298.15), P(1e5), charge(0),
      b(arma::zeros(system.elements().size()))
    {}
};

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumProblem::EquilibriumProblem(const EquilibriumProblem& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumProblem::~EquilibriumProblem()
{}

auto EquilibriumProblem::operator=(EquilibriumProblem other) -> EquilibriumProblem&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumProblem::setTemperature(double val) -> EquilibriumProblem&
{
    pimpl->T = val;
    return *this;
}

auto EquilibriumProblem::setPressure(double val) -> EquilibriumProblem&
{
    pimpl->P = val;
    return *this;
}

auto EquilibriumProblem::setElementAmounts(const Vector& b) -> EquilibriumProblem&
{
    pimpl->b = b;
    return *this;
}

auto EquilibriumProblem::setCharge(double val) -> EquilibriumProblem&
{
    pimpl->charge = val;
    return *this;
}

auto EquilibriumProblem::temperature() const -> double
{
    return pimpl->T;
}

auto EquilibriumProblem::pressure() const -> double
{
    return pimpl->P;
}

auto EquilibriumProblem::charge() const -> double
{
    return pimpl->charge;
}

auto EquilibriumProblem::elementAmounts() const -> const Vector&
{
    return pimpl->b;
}

auto EquilibriumProblem::balanceMatrix() const -> const Matrix&
{
    return pimpl->A;
}

auto EquilibriumProblem::independentComponents() const -> const Indices&
{
    return pimpl->independent_components;
}

auto EquilibriumProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

/// Create the objective function for the Gibbs energy minimization problem
auto createObjectiveFunction(const EquilibriumProblem& problem) -> ObjectiveFunction
{
    const double T = problem.temperature();
    const double P = problem.pressure();
    ChemicalVector u;
    ObjectiveResult res;

    ObjectiveFunction fn = [=](const Vector& n) mutable
    {
        u = problem.system().chemicalPotentials(T, P, n);
        res.func = arma::dot(n, u.val());
        res.grad = u.val();
        res.hessian = u.ddn();
        return res;
    };

    return fn;
}

/// Create constraint function for the Gibbs energy minimization problem
auto createConstraintFunction(const EquilibriumProblem& problem) -> ConstraintFunction
{
    // The right-hand side vector of the balance constraint
    Vector b = problem.elementAmounts();
    Vector z = arma::ones(1) * problem.charge();
    b = arma::join_vert(b, z);
    b = rows(b, problem.independentComponents());

    // The result of the equilibrium constraint evaluation
    ConstraintResult res;

    // Set the gradient of the constraint function as the formula matrix with linearly independent rows
    res.grad = problem.balanceMatrix();

    // Define the component (mass and charge) balance contraints
    ConstraintFunction fn = [=](const Vector& n) mutable
    {
        res.func = res.grad*n - b;
        return res;
    };

    return fn;
}

EquilibriumProblem::operator OptimumProblem() const
{
    const unsigned num_equilibrium_species = partition().equilibriumSpeciesIndices().size();
    const unsigned num_components = independentComponents().size();

    OptimumProblem problem(num_equilibrium_species, num_components);
    problem.setObjective(createObjectiveFunction(*this));
    problem.setConstraint(createConstraintFunction(*this));
    problem.setLowerBounds(0.0);

    return problem;
}

} // namespace Reaktor
