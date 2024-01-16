// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Common/MolalityUtils.hpp>
#include <Reaktoro/Common/MoleFractionUtils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

using autodiff::jacobian;
using autodiff::wrt;
using autodiff::at;

struct EquilibriumHessian::Impl
{
    /// The chemical system associated to this object.
    const ChemicalSystem system;

    /// The chemical properties of the system.
    ChemicalProps props;

    /// The auxiliary matrix for ∂(µ/RT)/∂n.
    MatrixXd dudn;

    /// The auxiliary vector to compute the diagonal of ∂(µ/RT)/∂n.
    VectorXd dudn_diag;

    /// The auxiliary vector of species amounts.
    VectorXr n;

    /// The functions for each phase that assemble the block of approximate derivatives in ∂(µ/RT)/∂n.
    Vec<Fn<void(VectorXrConstRef, MatrixXdRef)>> approxfuncs;

    /// The functions for each phase that assemble the diagonal of approximate derivatives in ∂(µ/RT)/∂n.
    Vec<Fn<void(VectorXrConstRef, VectorXdRef)>> approxfuncsdiag;

    ///
    Impl(ChemicalSystem const& system)
    : system(system), props(system)
    {
        const auto numphases = system.phases().size();
        const auto numspecies = system.species().size();

        dudn.resize(numspecies, numspecies);
        dudn_diag.resize(numspecies);

        approxfuncs.resize(numphases);
        approxfuncsdiag.resize(numphases);

        auto iphase = 0;

        for(auto iphase = 0; iphase < numphases; ++iphase)
        {
            auto const& phase = system.phase(iphase);

            if(phase.aggregateState() == AggregateState::Aqueous)
            {
                // IDEAL ACTIVITY OF WATER
                // \ln a_{w}=-\dfrac{1-x_{w}}{x_{w}}=-\dfrac{n_{\Sigma}-n_{w}}{n_{w}}
                // \dfrac{\partial\ln a_{w}}{\partial n_{i}}=\begin{cases}-\dfrac{1}{n_{w}} & i\neq w\\\dfrac{n_{\Sigma}-n_{w}}{n_{w}^{2}} & i=w\end{cases}
                const auto iH2O = phase.species().indexWithFormula("H2O");
                approxfuncs[iphase] = [=](VectorXrConstRef const& np, MatrixXdRef block)
                {
                    lnMolalitiesJacobian(np, iH2O, block);
                    const auto nH2O = np[iH2O];
                    const auto nsum = np.sum();
                    block.row(iH2O).array() = -1.0/nH2O;
                    block(iH2O, iH2O) = (nsum - nH2O)/(nH2O*nH2O);
                };
                approxfuncsdiag[iphase] = [=](VectorXrConstRef const& np, VectorXdRef segment)
                {
                    lnMolalitiesJacobianDiagonal(np, iH2O, segment);
                    const auto nH2O = np[iH2O];
                    const auto nsum = np.sum();
                    segment[iH2O] = (nsum - nH2O)/(nH2O*nH2O);
                };
            }
            else
            {
                approxfuncs[iphase] = [=](VectorXrConstRef const& np, MatrixXdRef block)
                {
                    lnMoleFractionsJacobian(np, block);
                };
                approxfuncsdiag[iphase] = [=](VectorXrConstRef const& np, VectorXdRef segment)
                {
                    lnMoleFractionsJacobianDiagonal(np, segment);
                };
            }
        }
    }

    auto exact(real const& T, real const& P, VectorXrConstRef const& nconst) -> MatrixXdConstRef
    {
        n = nconst;
        auto fn = [&](VectorXrConstRef const& n) -> VectorXr
        {
            props.update(T, P, n);
            return props.speciesChemicalPotentials();
        };
        const double RT = universalGasConstant * T;
        dudn.noalias() = jacobian(fn, wrt(n), at(n))/RT;
        return dudn;
    }

    auto partiallyExact(real const& T, real const& P, VectorXrConstRef const& nconst, VectorXlConstRef const& idxs) -> MatrixXdConstRef
    {
        n = nconst;
        auto fn = [&](VectorXrConstRef const& n) -> VectorXr
        {
            props.update(T, P, n);
            return props.speciesChemicalPotentials();
        };
        const double RT = universalGasConstant * T;
        dudn = approximate(n);
        dudn(Eigen::all, idxs) = jacobian(fn, wrt(n(idxs)), at(n))/RT;
        return dudn;
    }

    auto approximate(VectorXrConstRef const& n) -> MatrixXdConstRef
    {
        dudn.fill(0.0); // clear previous state of dudn
        const auto numphases = system.phases().size();
        auto offset = 0;
        for(auto i = 0; i < numphases; ++i)
        {
            const auto length = system.phase(i).species().size();
            const auto np = n.segment(offset, length);
            auto dupdnp = dudn.block(offset, offset, length, length);
            approxfuncs[i](np, dupdnp);
            offset += length;
        }
        return dudn;
    }

    auto diagonal(VectorXrConstRef const& n) -> MatrixXdConstRef
    {
        dudn.fill(0.0); // clear previous state of dudn
        const auto numphases = system.phases().size();
        auto offset = 0;
        for(auto i = 0; i < numphases; ++i)
        {
            const auto length = system.phase(i).species().size();
            const auto np = n.segment(offset, length);
            const auto dupdnp_diag = dudn_diag.segment(offset, length);
            approxfuncsdiag[i](np, dupdnp_diag);
            offset += length;
        }
        dudn.diagonal() = dudn_diag;
        return dudn;
    }
};

EquilibriumHessian::EquilibriumHessian(ChemicalSystem const& system)
: pimpl(new Impl(system))
{
}

EquilibriumHessian::EquilibriumHessian(EquilibriumHessian const& other)
: pimpl(new Impl(*other.pimpl))
{
}

EquilibriumHessian::~EquilibriumHessian()
{
}

auto EquilibriumHessian::operator=(EquilibriumHessian other) -> EquilibriumHessian&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumHessian::exact(real const& T, real const& P, VectorXrConstRef const& n) -> MatrixXdConstRef
{
    return pimpl->exact(T, P, n);
}

auto EquilibriumHessian::partiallyExact(real const& T, real const& P, VectorXrConstRef const& n, VectorXlConstRef const& idxs) -> MatrixXdConstRef
{
    return pimpl->partiallyExact(T, P, n, idxs);
}

auto EquilibriumHessian::approximate(VectorXrConstRef const& n) -> MatrixXdConstRef
{
    return pimpl->approximate(n);
}

auto EquilibriumHessian::diagonal(VectorXrConstRef const& n) -> MatrixXdConstRef
{
    return pimpl->diagonal(n);
}

} // namespace Reaktoro
