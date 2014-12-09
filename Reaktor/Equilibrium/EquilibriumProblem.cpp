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

    /// The right-hand side vector of the balance equations
    Matrix bz;

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

auto EquilibriumProblem::constraint() const -> EquilibriumConstraint
{
    Vector b = elementAmounts();
    Vector z = arma::ones(1) * charge();
    b = arma::join_rows(b, z);
    b = rows(b, independentComponents());

    EquilibriumConstraintResult res;
    res.ddn = balanceMatrix();

    auto func = [=](const Vector& n) mutable
    {
        res.val = res.ddn*n - b;
        return res;
    };

    return func;
}

} // namespace Reaktor
