// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "EquilibriumJacobian.hpp"

// Reaktoro includes
#include <Reaktoro/Common/MolalityUtils.hpp>
#include <Reaktoro/Common/MoleFractionUtils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

using autodiff::jacobian;
using autodiff::wrt;
using autodiff::at;

struct EquilibriumJacobian::Impl
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
    Impl(const ChemicalSystem& system)
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
            const auto& phase = system.phase(iphase);

            if(phase.aggregateState() == AggregateState::Aqueous)
            {
                const auto iH2O = phase.species().indexWithFormula("H2O");
                approxfuncs[iphase] = [=](VectorXrConstRef np, MatrixXdRef block)
                {
                    lnMolalitiesJacobian(np, iH2O, block);
                };
                approxfuncsdiag[iphase] = [=](VectorXrConstRef np, VectorXdRef segment)
                {
                    lnMolalitiesJacobianDiagonal(np, iH2O, segment);
                };
            }
            else
            {
                approxfuncs[iphase] = [=](VectorXrConstRef np, MatrixXdRef block)
                {
                    lnMoleFractionsJacobian(np, block);
                };
                approxfuncsdiag[iphase] = [=](VectorXrConstRef np, VectorXdRef segment)
                {
                    lnMoleFractionsJacobianDiagonal(np, segment);
                };
            }
        }
    }

    auto dudnExact(const real& T, const real& P, VectorXrConstRef nconst) -> MatrixXdConstRef
    {
        n = nconst;
        auto fn = [&](VectorXrConstRef n) -> VectorXr
        {
            props.update(T, P, n);
            return props.speciesChemicalPotentials();
        };
        const double RT = universalGasConstant * T;
        dudn.noalias() = jacobian(fn, wrt(n), at(n))/RT;
        return dudn;
    }

    auto dudnPartiallyExact(const real& T, const real& P, VectorXrConstRef nconst, VectorXlConstRef idxs) -> MatrixXdConstRef
    {
        n = nconst;
        auto fn = [&](VectorXrConstRef n) -> VectorXr
        {
            props.update(T, P, n);
            return props.speciesChemicalPotentials();
        };
        const double RT = universalGasConstant * T;
        dudn = dudnApproximate(n);
        dudn(Eigen::all, idxs) = jacobian(fn, wrt(n(idxs)), at(n))/RT;
        return dudn;
    }

    auto dudnApproximate(VectorXrConstRef n) -> MatrixXdConstRef
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

    auto dudnDiagonal(VectorXrConstRef n) -> MatrixXdConstRef
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

EquilibriumJacobian::EquilibriumJacobian(const ChemicalSystem& system)
: pimpl(new Impl(system))
{
}

EquilibriumJacobian::EquilibriumJacobian(const EquilibriumJacobian& other)
: pimpl(new Impl(*other.pimpl))
{
}

EquilibriumJacobian::~EquilibriumJacobian()
{
}

auto EquilibriumJacobian::operator=(EquilibriumJacobian other) -> EquilibriumJacobian&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumJacobian::dudnExact(const real& T, const real& P, VectorXrConstRef n) -> MatrixXdConstRef
{
    return pimpl->dudnExact(T, P, n);
}

auto EquilibriumJacobian::dudnPartiallyExact(const real& T, const real& P, VectorXrConstRef n, VectorXlConstRef idxs) -> MatrixXdConstRef
{
    return pimpl->dudnPartiallyExact(T, P, n, idxs);
}

auto EquilibriumJacobian::dudnApproximate(VectorXrConstRef n) -> MatrixXdConstRef
{
    return pimpl->dudnApproximate(n);
}

auto EquilibriumJacobian::dudnDiagonal(VectorXrConstRef n) -> MatrixXdConstRef
{
    return pimpl->dudnDiagonal(n);
}

} // namespace Reaktoro
