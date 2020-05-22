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
class ChemicalProps;
class ChemicalSystem;
class EquilibriumConstraints;

/// The objective function to be minimized in a chemical equilibrium calculation.
struct EquilibriumObjective
{
    /// The function that computes the value of the objective function.
    Fn<real(const ChemicalProps&)> f;

    /// The function that computes the gradient vector of the objective function.
    Fn<void(const ChemicalProps&, VectorXrRef)> g;

    /// The function that computes the Hessian matrix of the objective function.
    Fn<void(const ChemicalProps&, MatrixXdRef)> H;
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

    /// Return the number of components associated with given equilibrium constraints.
    /// The number of components is the sum of the number of elements in the
    /// system, including electric charge, and the number of inert reactions.
    auto numComponents() const -> Index;

    /// Return the total number of variables associated with given equilibrium constraints.
    /// The number of variables is the sum of the number of species, the number
    /// of chemical potential constraints, and the number of introduced control
    /// variables such as temperature, pressure and/or the amounts of titrants.
    auto numVariables() const -> Index;

    /// Return the total number of introduced control variables associated with given equilibrium constraints.
    /// The control variables are those, which must match with the number of imposed functional constraints, is the number of sum of the number of species, the number
    /// of chemical potential constraints, and the number of introduced control
    /// variables such as temperature, pressure and/or the amounts of titrants.
    auto numControlVariables() const -> Index;

    /// Assemble the conservation matrix based on the given equilibrium constraints.
    auto conservationMatrix() const -> MatrixXd;

    /// Assemble the objective function to be minimized based on the given equilibrium constraints.
    auto objective() const -> EquilibriumObjective;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Write here a brief about this class
class GibbsHessian
{
public:
    /// Construct a default GibbsHessian object.
    GibbsHessian();

private:

};


} // namespace Reaktoro
