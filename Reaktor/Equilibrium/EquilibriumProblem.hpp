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

    /// Set the molar amounts of each element for the equilibrium calculation (in units of mol)
    /// @param b The vector of molar amounts of each element (in units of mol)
    auto setElementAmounts(const Vector& b) -> EquilibriumProblem&;

    /// Set the molar amounts of each element for the equilibrium calculation (in units of mol)
    /// @param amount The same molar amount for the all elements (in units of mol)
    auto setElementAmounts(double amount) -> EquilibriumProblem&;

    /// Add a given amount of a compound or species to the equilibrium recipe.
    /// This method will first check if the given compound is present in the chemical system.
    /// If true, then `addSpecies` will be called. Otherwise, `addCompound` will be called.
    /// @param name The name of the compound or species
    /// @param amount The amount of the compound or species
    /// @param units The units of the amount (must be convertible to either mol or kg)
    auto add(std::string name, double amount, std::string units) -> EquilibriumProblem&;

    /// Add a given amount of a compound to the equilibrium recipe.
    /// The compound must not have a chemical element that is not present in the chemical system.
    /// @param name The name of the compound (e.g., H2O, CaCO3)
    /// @param amount The amount of the compound
    /// @param units The units of the amount (must be convertible to either mol or kg)
    auto addCompound(std::string name, double amount, std::string unit) -> EquilibriumProblem&;

    /// Add a given amount of a species to the equilibrium recipe
    /// The species must be present in the chemical system.
    /// @param name The name of the species
    /// @param amount The amount of the species
    /// @param units The units of the amount (must be convertible to either mol or kg)
    auto addSpecies(std::string name, double amount, std::string unit) -> EquilibriumProblem&;

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

    /// The balance matrix of the chemical system
    auto balanceMatrix() const -> const Matrix&;

    /// Get a reference to the ChemicalSystem instance used to create this EquilibriumProblem instance
    auto system() const -> const ChemicalSystem&;

    /// Get a reference to the Partition instance used to create this EquilibriumProblem instance
    auto partition() const -> const Partition&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
