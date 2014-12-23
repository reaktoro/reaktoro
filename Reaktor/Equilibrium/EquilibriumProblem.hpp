// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalSystem;
class Partition;
class OptimumProblem;

/// A type that defines an equilibrium problem
class EquilibriumProblem
{
public:
    /// Construct an EquilibriumProblem instance
    explicit EquilibriumProblem(const ChemicalSystem& system);

    /// Construct an EquilibriumProblem instance with partitioning information
    EquilibriumProblem(const ChemicalSystem& system, const Partition& partition);

    /// Construct a copy of a EquilibriumProblem instance
    EquilibriumProblem(const EquilibriumProblem& other);

    /// Destroy this EquilibriumProblem instance
    virtual ~EquilibriumProblem();

    /// Assign an EquilibriumProblem instance to this
    auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

    /// Set the temperature for the equilibrium calculation (in units of K)
    /// By default, the temperature is 298.15 K.
    auto setTemperature(double val) -> EquilibriumProblem&;

    /// Set the pressure for the equilibrium calculation (in units of Pa)
    /// By default, the pressure is 10<sup>5</sup> Pa.
    auto setPressure(double val) -> EquilibriumProblem&;

    /// Set the electrical charge for the equilibrium calculation (in units of mol)
    /// By default, the electrical charge is zero.
    auto setCharge(double val) -> EquilibriumProblem&;

    /// Set the amounts of the elements for the equilibrium calculation (in units of mol)
    auto setElementAmounts(const Vector& b) -> EquilibriumProblem&;

    /// Get the temperature for the equilibrium calculation (in units of K)
    auto temperature() const -> double;

    /// Get the pressure for the equilibrium calculation (in units of Pa)
    auto pressure() const -> double;

    /// Get the electrical charge for the equilibrium calculation (in units of mol)
    auto charge() const -> double;

    /// Get the amounts of the elements for the equilibrium calculation (in units of mol)
    auto elementAmounts() const -> const Vector&;

    /// Get the amounts of the components (elements and charge) for the equilibrium calculation (in units of mol)
    auto componentAmounts() const -> Vector;

    /// The balance matrix of the chemical system with linearly independent rows
    auto balanceMatrix() const -> const Matrix&;

    /// The indices of the linearly independent components
    auto components() const -> const Indices&;

    /// Get a reference to the ChemicalSystem instance used to create this EquilibriumProblem instance
    auto system() const -> const ChemicalSystem&;

    /// Get a reference to the Partition instance used to create this EquilibriumProblem instance
    auto partition() const -> const Partition&;

    /// Convert this EquilibriumProblem instance into a OptimumProblem instance
    operator OptimumProblem() const;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
