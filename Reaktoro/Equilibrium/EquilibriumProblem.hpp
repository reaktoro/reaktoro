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
class EquilibriumConstraints;
class EquilibriumDims;

/// The objective function to be minimized in a chemical equilibrium calculation.
struct EquilibriumObjective
{
    /// The function that computes the value of the objective function.
    Fn<double(VectorXdConstRef x)> f;

    /// The function that computes the gradient vector of the objective function.
    Fn<void(VectorXdConstRef x, VectorXdRef res)> g;

    /// The function that computes the Hessian matrix of the objective function.
    Fn<void(VectorXdConstRef x, MatrixXdRef res)> H;
};

/// The class used to define an equilibrium problem.
class EquilibriumProblem
{
public:
    /// Construct an EquilibriumProblem object with given constraints.
    explicit EquilibriumProblem(const EquilibriumConstraints& constraints);

    /// Construct a copy of an EquilibriumProblem object.
    EquilibriumProblem(const EquilibriumProblem& other);

    /// Destroy this EquilibriumProblem object.
    ~EquilibriumProblem();

    /// Assign a copy of an EquilibriumProblem object to this.
    auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

    /// Update the equilibrium constraints for the next chemical equilibrium calculation.
    /// @warning An error will result if new constraints are imposed. This
    /// method should only be used to update the parameters in the existing
    /// constraints. For example, a calculation was done before with a fixed pH
    /// value, and a new value needs to be imposed for the next calculation.
    auto update(const EquilibriumConstraints& constraints) -> void;

    /// Return the dimensions of the variables in the equilibrium problem.
    auto dims() const -> const EquilibriumDims&;

    /// Return the conservation matrix based on the given equilibrium constraints.
    auto conservationMatrix() const -> MatrixXd;

    /// Return the objective function to be minimized based on the given equilibrium constraints.
    /// @param state0 The initial chemical state of the system.
    auto objective(const ChemicalState& state0) const -> EquilibriumObjective;

    /// Set the lower bounds of the species amounts based on the given equilibrium constraints.
    /// @param state0 The initial chemical state of the system.
    /// @param res The array where the lower bounds are set.
    auto xlower(const ChemicalState& state0, ArrayXdRef res) const -> void;

    /// Set the upper bounds of the species amounts based on the given equilibrium constraints.
    /// @param state0 The initial chemical state of the system.
    /// @param res The array where the upper bounds are set.
    auto xupper(const ChemicalState& state0, ArrayXdRef res) const -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
