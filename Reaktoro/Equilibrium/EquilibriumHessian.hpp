// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
class ChemicalProps;
class ChemicalSystem;

/// The class used to compute the Hessian matrix of the Gibbs energy function.
class EquilibriumHessian
{
public:
    /// Construct an EquilibriumHessian object.
    EquilibriumHessian(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumHessian object.
    EquilibriumHessian(const EquilibriumHessian& other);

    /// Destroy this EquilibriumHessian object.
    ~EquilibriumHessian();

    /// Assign a copy of an EquilibriumHessian object to this.
    auto operator=(EquilibriumHessian other) -> EquilibriumHessian&;

    /// Return the exact Hessian matrix of the Gibbs energy function.
    /// This method computes exact molar derivatives of the chemical potentials
    /// for all species in the chemical system.
    auto exact(const ChemicalProps& props) -> MatrixXdConstRef;

    /// Return a partially exact Hessian matrix of the Gibbs energy function.
    /// This method computes exact molar derivatives of the chemical potentials
    /// for those species with given `indices` only. For all other species, and
    /// approximation is made based on ideal activities, essentially, using
    /// molar derivatives of mole fractions.
    /// @note This is the appropriate method for chemical equilibrium
    /// calculations in which several or many species have very low amounts and
    /// play no major role in the chemical equilibrium state being calculated.
    /// Approximate derivatives for these species suffice for numerical
    /// convergence reasons.
    auto partiallyExact(const ChemicalProps& props, const Indices& indices) -> MatrixXdConstRef;

    /// Return an approximation of the Hessian matrix of the Gibbs energy function.
    /// This method computes an approximation for the Hessian matrix based on
    /// ideal activities, essentially, using molar derivatives of mole
    /// fractions. This method is suitable when the chemical equilibrium
    /// calculation is converging towards a solution without difficulties.
    auto approx(const ChemicalProps& props) -> MatrixXdConstRef;

    /// Return a diagonal approximation of the Hessian matrix of the Gibbs energy function.
    /// This method computes a diagonal approximation for the Hessian matrix
    /// based on ideal activities, essentially, using molar derivatives of mole
    /// fractions. This method is suitable when the chemical equilibrium
    /// calculation is converging towards a solution without difficulties,
    /// resulting in very efficient linear algebra in each Newton stepping.
    auto diagonalApprox(const ChemicalProps& props) -> MatrixXdConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
