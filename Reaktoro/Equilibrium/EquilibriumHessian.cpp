// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "EquilibriumHessian.hpp"

// Reaktoro includes
#include <Reaktoro/Common/AutoDiff.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/MoleFractionUtils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

struct EquilibriumHessian::Impl
{
    /// The chemical system for which Hessian matrix is computed.
    ChemicalSystem system;

    /// The chemical properties of the system at which the Hessian needs to be computed.
    ChemicalProps props;

    /// The Hessian matrix.
    MatrixXd H;

    /// The number of species in the system.
    Index N;

    /// The temperature of the system (in K).
    real T = 0.0;

    /// The pressure of the system (in Pa).
    real P = 0.0;

    /// The amounts of each species in the system (in mol).
    ArrayXr n;

    /// Construct an EquilibriumHessian::Impl object.
    Impl(const ChemicalSystem& system)
    : system(system), props(system), N(system.species().size())
    {
        H = MatrixXd::Zero(N, N);
    }

    /// Update the temperature, pressure, and species amounts.
    auto update(const ChemicalProps& props0)
    {
        T = props.temperature();
        P = props.pressure();
        n = props.speciesAmounts();
    }

    /// Return the exact Hessian matrix of the Gibbs energy function.
    auto exact(const ChemicalProps& props0) -> MatrixXdConstRef
    {
        update(props);
        for(auto i = 0; i < N; ++i) {
            autodiff::seed(n[i]);
            props.update(T, P, n);
            autodiff::unseed(n[i]);
            H.row(i) = grad(props.chemicalPotentials());
        }
        return H;
    }

    /// Return a partially exact Hessian matrix of the Gibbs energy function.
    auto partiallyExact(const ChemicalProps& props0, const Indices& indices) -> MatrixXdConstRef
    {
        /// Compute the first layer of H with approximate derivatives.
        approx(props);

        /// Overwrite the approximations for those species with given indices.
        update(props);
        for(auto i : indices)
        {
            autodiff::seed(n[i]);
            props.update(T, P, n);
            autodiff::unseed(n[i]);
            H.row(i) = grad(props.chemicalPotentials());
        }

        return H;
    }

    /// Return an approximation of the Hessian matrix of the Gibbs energy function.
    auto approx(const ChemicalProps& props0) -> MatrixXdConstRef
    {
        const auto T = props0.temperature();
        const auto n = props0.speciesAmounts();
        const auto R = universalGasConstant;

        const double RT = R*T; // the use of double is intentional!

        H.fill(0.0); // clear previous state of H

        const auto Np = system.phases().size();
        auto offset = 0;
        for(auto i = 0; i < Np; ++i)
        {
            const auto length = system.phase(i).species().size();
            const auto np = n.segment(offset, length);
            MatrixXd Hp = H.block(offset, offset, length, length);
            lnMoleFractionsJacobian(np, Hp); // set Hp = d(lnx)/dn
            Hp *= RT; // finally set Hp = RT*d(lnx)/dn
            offset += length;
        }

        return H;
    }

    /// Return a diagonal approximation of the Hessian matrix of the Gibbs energy function.
    auto diagonalApprox(const ChemicalProps& props0) -> MatrixXdConstRef
    {
        const auto T = props0.temperature();
        const auto x = props0.moleFractions();
        const auto n = props0.speciesAmounts();
        const auto R = universalGasConstant;

        const double RT = R*T; // the use of double is intentional!

        H.fill(0.0); // clear previous state of H

        for(auto i = 0; i < N; ++i)
            H(i, i) = RT * (1 - x[i])/n[i];

        return H;
    }
};

EquilibriumHessian::EquilibriumHessian(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumHessian::EquilibriumHessian(const EquilibriumHessian& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumHessian::~EquilibriumHessian()
{}

auto EquilibriumHessian::operator=(EquilibriumHessian other) -> EquilibriumHessian&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumHessian::exact(const ChemicalProps& props) -> MatrixXdConstRef
{
    return pimpl->exact(props);
}

auto EquilibriumHessian::partiallyExact(const ChemicalProps& props, const Indices& indices) -> MatrixXdConstRef
{
    return pimpl->partiallyExact(props, indices);
}

auto EquilibriumHessian::approx(const ChemicalProps& props) -> MatrixXdConstRef
{
    return pimpl->approx(props);
}

auto EquilibriumHessian::diagonalApprox(const ChemicalProps& props) -> MatrixXdConstRef
{
    return pimpl->diagonalApprox(props);
}

} // namespace Reaktoro
