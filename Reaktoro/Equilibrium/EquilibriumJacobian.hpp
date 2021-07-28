// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// Used to compute the Jacobian of chemical potentials during the equilibrium calculations.
class EquilibriumJacobian
{
public:
    /// Construct an EquilibriumJacobian object with given specifications.
    explicit EquilibriumJacobian(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumJacobian object.
    EquilibriumJacobian(const EquilibriumJacobian& other);

    /// Destroy this EquilibriumJacobian object.
    ~EquilibriumJacobian();

    /// Assign a copy of an EquilibriumJacobian object to this.
    auto operator=(EquilibriumJacobian other) -> EquilibriumJacobian&;

    /// Evaluate the Jacobian matrix *∂(µ/RT)/∂n* with exact derivatives. This
    /// method uses automatic differentiation to compute the exact derivatives
    /// of *µ/RT* with respect to all species in the system.
    auto dudnExact(const real& T, const real& P, VectorXrConstRef n) -> MatrixXdConstRef;

    /// Evaluate the Jacobian matrix *∂(µ/RT)/∂n* with exact derivatives for
    /// selected species. This method uses automatic differentiation to compute
    /// the exact derivatives of *µ/RT* with respect to the species whose indices
    /// are given in `idxs`. For all other species, approximate derivatives are
    /// used instead. These approximate derivatives are obtained from ideal
    /// thermodynamic models for the phases. Thus, with this method, the
    /// derivatives in *∂(µ/RT)/∂n* associated with species with low amounts,
    /// which can be a large number, are quickly evaluated.
    auto dudnPartiallyExact(const real& T, const real& P, VectorXrConstRef n, VectorXlConstRef idxs) -> MatrixXdConstRef;

    /// Evaluate the Jacobian matrix *∂(µ/RT)/∂n* with approximate derivatives.
    /// These derivatives are analytically (and thus quickly) computed from
    /// ideal thermodynamic models for the phases in the system. For example,
    /// for liquid, gaseous, and solid solutions *∂(µ/RT)/∂n ≡ ∂(ln(x))/∂n*,
    /// where *x* are species mole fractions. For an aqueous solution, however,
    /// *∂(µ/RT)/∂n ≡ ∂(ln(m))/∂n*, where *m* are species molalities.
    auto dudnApproximate(VectorXrConstRef n) -> MatrixXdConstRef;

    /// Evaluate the Jacobian matrix *∂(µ/RT)/∂n* as a diagonal matrix using
    /// approximate derivatives. The computed diagonal matrix with this
    /// function is equivalent to extracting the diagonal entries from the
    /// matrix produced with @ref dudnApproximate.
    auto dudnDiagonal(VectorXrConstRef n) -> MatrixXdConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
