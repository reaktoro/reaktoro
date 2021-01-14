// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #pragma once

// // Optima includes
// #include <Optima/Problem.hpp>
// #include <Optima/State.hpp>

// // Reaktoro includes
// #include <Reaktoro/Common/Types.hpp>
// #include <Reaktoro/Common/Matrix.hpp>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalState;
// class ChemicalSystem;
// class EquilibriumSpecs;
// class EquilibriumDims;
// class EquilibriumOptions;

// /// The objective function to be minimized in a chemical equilibrium calculation.
// struct EquilibriumObjective
// {
//     /// The function that computes the value of the objective function.
//     Fn<double(VectorXdConstRef x)> f;

//     /// The function that computes the gradient vector of the objective function.
//     Fn<void(VectorXdConstRef x, VectorXdRef res)> g;

//     /// The function that computes the Hessian matrix of the objective function.
//     Fn<void(VectorXdConstRef x, MatrixXdRef res)> H;
// };

// /// The class used to define an equilibrium problem.
// class EquilibriumProblem
// {
// public:
//     /// Construct an EquilibriumProblem object with given specifications.
//     explicit EquilibriumProblem(const EquilibriumSpecs& specs);

//     /// Construct a copy of an EquilibriumProblem object.
//     EquilibriumProblem(const EquilibriumProblem& other);

//     /// Destroy this EquilibriumProblem object.
//     ~EquilibriumProblem();

//     /// Assign a copy of an EquilibriumProblem object to this.
//     auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

//     /// Set the options for the solution of the equilibrium problem.
//     auto setOptions(const EquilibriumOptions& options) -> void;

//     /// Return the dimensions of the variables in the equilibrium problem.
//     auto dims() const -> const EquilibriumDims&;

//     // /// Return the objective function to be minimized based on the given equilibrium conditions.
//     // /// @param state0 The initial chemical state of the system.
//     // auto objective(const ChemicalState& state0) const -> EquilibriumObjective;

//     // /// Set the lower bounds of the species amounts based on the given equilibrium conditions.
//     // /// @param state0 The initial chemical state of the system.
//     // /// @param res The array where the lower bounds are set.
//     // auto xlower(const ChemicalState& state0, ArrayXdRef res) const -> void;

//     // /// Set the upper bounds of the species amounts based on the given equilibrium conditions.
//     // /// @param state0 The initial chemical state of the system.
//     // /// @param res The array where the upper bounds are set.
//     // auto xupper(const ChemicalState& state0, ArrayXdRef res) const -> void;

//     /// Assemble the matrix @eq{A_{\mathrm{ex}}} in optimization problem.
//     auto assembleMatrixAex() const -> MatrixXd;

//     /// Assemble the matrix @eq{A_{\mathrm{ep}}} in optimization problem.
//     auto assembleMatrixAep() const -> MatrixXd;

//     /// Evaluate the objective function of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalObjectiveValue(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> real;

//     /// Evaluate the gradient of the objective function of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalObjectiveGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef;

//     /// Evaluate the Hessian of the objective function of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalObjectiveHessianX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

//     /// Evaluate the Hessian of the objective function of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalObjectiveHessianP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

//     /// Evaluate the equation constraints of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalEquationConstraints(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> VectorXrConstRef;

//     /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalEquationConstraintsGradX(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

//     /// Evaluate the Jacobian of the equation constraints of the chemical equilibrium problem.
//     /// @param x The amounts of the species and implicit titrants, @eq{x = (n, q)}.
//     /// @param p The amounts of the explicit titrants in the chemical equilibrium problem.
//     /// @param params The input parameters in the chemical equilibrium problem.
//     auto evalEquationConstraintsGradP(VectorXrConstRef x, VectorXrConstRef p, const Params& params) -> MatrixXdConstRef;

//     // auto createOptimaProblem() const -> Optima::Problem;

//     // auto createOptimaState() const -> Optima::State;

// private:
//     struct Impl;

//     Ptr<Impl> pimpl;
// };

// } // namespace Reaktoro
