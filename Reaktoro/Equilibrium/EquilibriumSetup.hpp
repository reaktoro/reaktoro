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
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumDims;
class EquilibriumProps;
class EquilibriumRestrictions;
class EquilibriumSpecs;
struct EquilibriumOptions;

/// The class used to construct the optimization problem for a chemical equilibrium calculation.
class EquilibriumSetup
{
public:
    /// Construct an EquilibriumSetup object with given specifications.
    explicit EquilibriumSetup(const EquilibriumSpecs& specs);

    /// Construct a copy of an EquilibriumSetup object.
    EquilibriumSetup(const EquilibriumSetup& other);

    /// Destroy this EquilibriumSetup object.
    ~EquilibriumSetup();

    /// Assign a copy of an EquilibriumSetup object to this.
    auto operator=(EquilibriumSetup other) -> EquilibriumSetup&;

    /// Set the options for the solution of the equilibrium problem.
    auto setOptions(const EquilibriumOptions& options) -> void;

    /// Return the dimensions of the variables in the equilibrium problem.
    auto dims() const -> const EquilibriumDims&;

    /// Return the options for the solution of the equilibrium problem.
    auto options() const -> const EquilibriumOptions&;

    /// Assemble the coefficient matrix `Aex` in the optimization problem.
    auto assembleMatrixAex() const -> MatrixXd;

    /// Assemble the coefficient matrix `Aep` in the optimization problem.
    auto assembleMatrixAep() const -> MatrixXd;

    /// Assemble the lower bound vector `xlower` in the optimization problem where *x = (n, q)*.
    /// @param restrictions The lower and upper bounds information of the species.
    /// @param state0 The initial chemical state of the system.
    auto assembleLowerBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd;

    /// Assemble the upper bound vector `xupper` in the optimization problem where *x = (n, q)*.
    /// @param restrictions The lower and upper bounds information of the species.
    /// @param state0 The initial chemical state of the system.
    auto assembleUpperBoundsVector(const EquilibriumRestrictions& restrictions, const ChemicalState& state0) const -> VectorXd;

    /// Evaluate the objective function of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> real;

    /// Evaluate the gradient of the objective function of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> VectorXrConstRef;

    /// Evaluate the Hessian of the objective function of the chemical equilibrium problem with respect to @p x and @p x.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalObjectiveHessianX(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Evaluate the Hessian of the objective function of the chemical equilibrium problem with respect to @p x and @p p.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Evaluate the Hessian of the objective function of the chemical equilibrium problem with respect to @p x and @p w.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalObjectiveHessianW(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Evaluate the equation constraints of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> VectorXrConstRef;

    /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem with respect to @p x.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalEquationConstraintsGradX(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem with respect to @p p.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem with respect to @p w.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param w The input variables *w* in the chemical equilibrium problem.
    auto evalEquationConstraintsGradW(VectorXrConstRef x, VectorXrConstRef p, VectorXrConstRef w) -> MatrixXdConstRef;

    /// Enable recording of derivatives of the chemical properties with respect
    /// to *(n, p, w)* to contruct its full Jacobian matrix.
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
    auto equilibriumProps() const -> const EquilibriumProps&;

    /// Return the current chemical properties of the system as a ChemicalProps object.
    auto chemicalProps() const -> const ChemicalProps&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
