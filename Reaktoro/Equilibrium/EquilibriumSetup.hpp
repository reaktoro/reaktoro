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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumDims;
class EquilibriumOptions;
class EquilibriumRestrictions;
class EquilibriumSpecs;
class Params;

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

    /// Assemble the right-hand side vector `be` in the optimization problem.
    auto assembleVectorBe(const EquilibriumConditions& conditions, const ChemicalState& state0) -> VectorXr;

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
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> real;

    /// Evaluate the gradient of the objective function of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef;

    /// Evaluate the Hessian of the objective function of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalObjectiveHessianX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

    /// Evaluate the Hessian of the objective function of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

    /// Evaluate the equation constraints of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef;

    /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalEquationConstraintsGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

    /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem.
    /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

    /// Extract temperature from either vector `p` or `params`.
    /// If temperature in unknown in the specifications of the equilibrium
    /// problem, then it is extracted from `p`. Otherwise, temperature is known
    /// and available in the Params object `params`.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto extractTemperature(VectorXrConstRef p, const Params& params) const -> real;

    /// Extract pressure from either vector `p` or `params`.
    /// If pressure in unknown in the specifications of the equilibrium
    /// problem, then it is extracted from `p`. Otherwise, pressure is known
    /// and available in the Params object `params`.
    /// @param p The values of the *p* control variables (e.g., temperature, pressure, and/or amounts of explicit titrants).
    /// @param params The input parameters in the chemical equilibrium problem.
    auto extractPressure(VectorXrConstRef p, const Params& params) const -> real;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
