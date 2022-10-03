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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// Used to compute the Hessian matrix of the Gibbs energy function.
class EquilibriumHessian
{
public:
    /// Construct an EquilibriumHessian object with given specifications.
    explicit EquilibriumHessian(ChemicalSystem const& system);

    /// Construct a copy of an EquilibriumHessian object.
    EquilibriumHessian(EquilibriumHessian const& other);

    /// Destroy this EquilibriumHessian object.
    ~EquilibriumHessian();

    /// Assign a copy of an EquilibriumHessian object to this.
    auto operator=(EquilibriumHessian other) -> EquilibriumHessian&;

    /// Evaluate the Hessian matrix *∂(µ/RT)/∂n* with exact derivatives. This
    /// method uses automatic differentiation to compute the exact derivatives
    /// of *µ/RT* with respect to all species in the system.
    auto exact(real const& T, real const& P, VectorXrConstRef const& n) -> MatrixXdConstRef;

    /// Evaluate the Hessian matrix *∂(µ/RT)/∂n* with exact derivatives for selected species. This
    /// method uses automatic differentiation to compute the exact derivatives of *µ/RT* with
    /// respect to the species whose indices are given in `idxs`. For all other species, approximate
    /// derivatives are used instead. These approximate derivatives are obtained from ideal
    /// thermodynamic models for the phases. Thus, with this method, the derivatives in *∂(µ/RT)/∂n*
    /// associated with species with low amounts, which can be a large number, are quickly
    /// evaluated. This is appropriate for chemical equilibrium calculations in which several or
    /// many species have very low amounts and play no major role in the chemical equilibrium state
    /// being calculated. Approximate derivatives for these species suffice for numerical convergence.
    auto partiallyExact(real const& T, real const& P, VectorXrConstRef const& n, VectorXlConstRef const& idxs) -> MatrixXdConstRef;

    /// Evaluate the Hessian matrix *∂(µ/RT)/∂n* with approximate derivatives. These derivatives
    /// are analytically (and thus quickly) computed from ideal thermodynamic models for the phases
    /// in the system. For example, for liquid, gaseous, and solid solutions *∂(µ/RT)/∂n ≡ ∂(ln(x))/∂n*,
    /// where *x* are species mole fractions. For an aqueous solution,
    /// however, *∂(µ/RT)/∂n ≡ ∂(ln(m))/∂n*, where *m* are species molalities.
    auto approximate(VectorXrConstRef const& n) -> MatrixXdConstRef;

    /// Evaluate the Hessian matrix *∂(µ/RT)/∂n* as a diagonal matrix using approximate
    /// derivatives. The computed diagonal matrix with this function is equivalent to extracting the
    /// diagonal entries from the matrix produced with @ref dudnApproximate.
    auto diagonal(VectorXrConstRef const& n) -> MatrixXdConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
