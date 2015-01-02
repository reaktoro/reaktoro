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

// C++ includes
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Constants.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>
#include <Reaktor/Math/MathUtils.hpp>
#include <Reaktor/Optimization/OptimumProblem.hpp>

namespace Reaktor {

struct EquilibriumProblem::Impl
{
    /// The reference to the ChemicalSystem instance
    ChemicalSystem system;

    /// The reference to the Partition instance
    Partition partition;

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
    Indices components;

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const ChemicalSystem& system)
    : Impl(system, Partition::allEquilibrium(system))
    {}

    /// Construct a EquilibriumProblem::Impl instance
    Impl(const ChemicalSystem& system, const Partition& partition)
    : system(system), partition(partition),
      T(298.15), P(1e5), charge(0.0), b(zeros(system.elements().size()))
    {
        // Set the balance matrix of the chemical system
        A = Reaktor::balanceMatrix(system);

        // Set the linearly independent components (the row indices of A that are linearly independent)
        components = linearlyIndependentRows(A, A);
    }
};

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumProblem::EquilibriumProblem(const ChemicalSystem& system, const Partition& partition)
: pimpl(new Impl(system, partition))
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

auto EquilibriumProblem::componentAmounts() const -> Vector
{
    const unsigned num_elements = system().elements().size();
    const unsigned num_components = components().size();
    Vector b(num_components);
    for(unsigned j = 0; j < num_components; ++j)
    {
        const Index icomponent = components()[j];
        b[j] = icomponent < num_elements ? elementAmounts()[icomponent] : charge();
    }
    return b;
}

auto EquilibriumProblem::balanceMatrix() const -> const Matrix&
{
    return pimpl->A;
}

auto EquilibriumProblem::components() const -> const Indices&
{
    return pimpl->components;
}

auto EquilibriumProblem::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto EquilibriumProblem::partition() const -> const Partition&
{
    return pimpl->partition;
}

EquilibriumProblem::operator OptimumProblem() const
{
    // The number of equilibrium species and linearly independent components in the equilibrium partition
    const unsigned num_equilibrium_species = partition().equilibriumSpeciesIndices().size();
    const unsigned num_components = components().size();

    // The reference to the chemical system
    const ChemicalSystem& system = pimpl->system;

    // The temperature and pressure of the equilibrium calculation
    const double T  = temperature();
    const double P  = pressure();
    const double RT = universalGasConstant*T;

    // An auxiliary vector to hold the chemical potentials of the species
    // scaled by RT. This is a shared pointer because (1) it must outlive
    // this function and (2) its content is shared among the functions below
    std::shared_ptr<Vector> u_ptr(new Vector());

    // Define the Gibbs energy function
    ObjectiveFunction gibbs = [=](const Vector& n) mutable
    {
        *u_ptr = (1/RT)*system.chemicalPotentials(T, P, n).val();
        return dot(n, *u_ptr);
    };

    // Define the gradient of the Gibbs energy function
    ObjectiveGradFunction gibbs_grad = [=](const Vector& n) mutable
    {
        return *u_ptr;
    };

    // Define the Hessian of the Gibbs energy function
    ObjectiveDiagonalHessianFunction gibbs_hessian = [=](const Vector& n) mutable
    {
        return inv(n);
    };

    // The left-hand side matrix of the linearly independent mass-charge balance equations
    const Matrix A = balanceMatrix();

    // The right-hand side vector of the linearly independent mass-charge balance equations
    const Vector b = componentAmounts();

    // Define the mass-cahrge balance contraint function
    ConstraintFunction balance_constraint = [=](const Vector& n) -> Vector
    {
        return A*n - b;
    };

    // Define the gradient function of the mass-cahrge balance contraint function
    ConstraintGradFunction balance_constraint_grad = [=](const Vector& n) -> Matrix
    {
        return A;
    };

    // Setup an OptimumProblem instance with the Gibbs energy function and the balance constraints
    OptimumProblem problem(num_equilibrium_species, num_components);
    problem.setObjective(gibbs);
    problem.setObjectiveGrad(gibbs_grad);
    problem.setObjectiveDiagonalHessian(gibbs_hessian);
    problem.setConstraint(balance_constraint);
    problem.setConstraintGrad(balance_constraint_grad);
    problem.setLowerBounds(0.0);

    return problem;
}

} // namespace Reaktor
