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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumDims;
class EquilibriumProps;
class EquilibriumRestrictions;
class EquilibriumSpecs;
struct EquilibriumOptions;

/// Used to construct the optimization problem for a chemical equilibrium calculation.
class EquilibriumSetup
{
public:
    /// Construct an EquilibriumSetup object with given specifications.
    explicit EquilibriumSetup(EquilibriumSpecs const& specs);

    /// Construct a copy of an EquilibriumSetup object.
    EquilibriumSetup(EquilibriumSetup const& other);

    /// Destroy this EquilibriumSetup object.
    ~EquilibriumSetup();

    /// Assign a copy of an EquilibriumSetup object to this.
    auto operator=(EquilibriumSetup other) -> EquilibriumSetup&;

    /// Set the options for the solution of the equilibrium problem.
    auto setOptions(EquilibriumOptions const& options) -> void;

    /// Return the dimensions of the variables in the equilibrium problem.
    auto dims() const -> EquilibriumDims const&;

    /// Return the options for the solution of the equilibrium problem.
    auto options() const -> EquilibriumOptions const&;

    /// Return the coefficient matrix `Aex` in the optimization problem.
    auto Aex() const -> MatrixXdConstRef;

    /// Return the coefficient matrix `Aep` in the optimization problem.
    auto Aep() const -> MatrixXdConstRef;

    /// Assemble the lower bound vector `xlower` in the optimization problem where *x = (n, q)*.
    /// @param restrictions The lower and upper bounds information of the species.
    /// @param state0 The initial chemical state of the system.
    auto assembleLowerBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd;

    /// Assemble the upper bound vector `xupper` in the optimization problem where *x = (n, q)*.
    /// @param restrictions The lower and upper bounds information of the species.
    /// @param state0 The initial chemical state of the system.
    auto assembleUpperBoundsVector(EquilibriumRestrictions const& restrictions, ChemicalState const& state0) const -> VectorXd;

    /// Update the chemical potentials and residuals of the equilibrium constraints.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto update(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> void;

    /// Update the derivatives of the chemical potentials and residuals of the equilibrium constraints with respect to *x*.
    /// @param ibasicvars The indices of the current basic variables in *x*.
    auto updateGradX(VectorXlConstRef ibasicvars) -> void;

    /// Update the derivatives of the chemical potentials and residuals of the equilibrium constraints with respect to *p*.
    auto updateGradP() -> void;

    /// Update the derivatives of the chemical potentials and residuals of the equilibrium constraints with respect to *w*.
    auto updateGradW() -> void;

    /// Get the updated Gibbs energy value.
    auto getGibbsEnergy() -> real;

    /// Get the updated gradient vector of the Gibbs energy with respect to *x*.
    auto getGibbsGradX() -> VectorXdConstRef;

    /// Get the updated Hessian matrix of the Gibbs energy with respect to *x*.
    auto getGibbsHessianX() -> MatrixXdConstRef;

    /// Get the updated Hessian matrix of the Gibbs energy with respect to *p*.
    auto getGibbsHessianP() -> MatrixXdConstRef;

    /// Get the updated Hessian matrix of the Gibbs energy with respect to *c*.
    auto getGibbsHessianC() -> MatrixXdConstRef;

    /// Get the updated residuals of the equilibrium constraints.
    auto getConstraintResiduals() -> MatrixXdConstRef;

    /// Get the updated Jacobian of the equilibrium constraints with respect to *x*.
    auto getConstraintResidualsGradX() -> MatrixXdConstRef;

    /// Get the updated Jacobian of the equilibrium constraints with respect to *p*.
    auto getConstraintResidualsGradP() -> MatrixXdConstRef;

    /// Get the updated Jacobian of the equilibrium constraints with respect to *c*.
    auto getConstraintResidualsGradC() -> MatrixXdConstRef;

    /// Return true if partially exact derivatives are adopted for the Hessian matrix *Hxx*.
    auto usingPartiallyExactDerivatives() -> bool;

    /// Return true if a diagonal structure is adopted for the Hessian matrix *Hxx*.
    auto usingDiagonalApproxDerivatives() -> bool;

    /// Enable recording of derivatives of the chemical properties with respect
    /// to *(n, p, w)* to construct its full Jacobian matrix.
    /// Consider a series of forward automatic differentiation passes to
    /// compute the partial derivatives of the chemical properties with respect
    /// to the variables *(n, p, w)*. Use this method before these operations
    /// so that these derivatives are recorded in this EquilibriumProps object.
    /// At the end of all forward passes, the full Jacobian of the chemical
    /// properties will have been constructed.
    /// @note Call @ref assembleChemicalPropsJacobianEnd after these forward passes
    /// have ended to eliminates the minor overhead of recording derivatives.
    auto assembleChemicalPropsJacobianBegin() -> void;

    /// Disable recording of derivatives of the chemical properties with
    /// respect to *(n, p, w)* to indicate the end of the full Jacobian matrix
    /// construction.
    auto assembleChemicalPropsJacobianEnd() -> void;

    /// Return the current chemical properties of the system as an EquilibriumProps object.
    auto equilibriumProps() const -> EquilibriumProps const&;

    /// Return the current chemical properties of the system as a ChemicalProps object.
    auto chemicalProps() const -> ChemicalProps const&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
