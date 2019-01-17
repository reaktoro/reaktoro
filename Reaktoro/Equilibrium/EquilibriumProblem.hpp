// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class Partition;
class Species;

/// A type that defines an equilibrium problem
class EquilibriumProblem
{
public:
    /// Construct an EquilibriumProblem instance
    explicit EquilibriumProblem(const ChemicalSystem& system);

    /// Construct a copy of a EquilibriumProblem instance
    EquilibriumProblem(const EquilibriumProblem& other);

    /// Destroy this EquilibriumProblem instance
    virtual ~EquilibriumProblem();

    /// Assign an EquilibriumProblem instance to this
    auto operator=(EquilibriumProblem other) -> EquilibriumProblem&;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(const Partition& partition) -> EquilibriumProblem&;

    /// Set the temperature for the equilibrium calculation.
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value (in units of K)
    auto setTemperature(double val) -> EquilibriumProblem&;

    /// Set the temperature for the equilibrium calculation with given units.
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value
    /// @param units The units of the temperature (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine)
    auto setTemperature(double val, std::string units) -> EquilibriumProblem&;

    /// Set the pressure for the equilibrium calculation.
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value (in units of Pa).
    auto setPressure(double val) -> EquilibriumProblem&;

    /// Set the pressure for the equilibrium calculation.
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value
    /// @param units The units of the pressure (Pa, kPa, MPa, GPa, atm, mmHg, inHg, psi, kpsi, Mpsi, psf, bar, torr, inH2O, ftH2O, pascal)
    auto setPressure(double val, std::string units) -> EquilibriumProblem&;

    /// Set the mole amounts of each element for the equilibrium calculation.
    /// @param b The vector of mole amounts of each element (in units of mol)
    auto setElementAmounts(VectorConstRef b) -> EquilibriumProblem&;

    /// Set the mole amounts of each element for the equilibrium calculation.
    /// @param amount The mole amount for all elements (in units of mol)
    auto setElementAmounts(double amount) -> EquilibriumProblem&;

    /// Set the mole amount of an element for the equilibrium calculation (in units of mol)
    /// @param ielement The index of the element
    /// @param amount The same mole amount for the all elements (in units of mol)
    auto setElementAmount(Index ielement, double amount) -> EquilibriumProblem&;

    /// Set the mole amount of an element for the equilibrium calculation (in units of mol)
    /// @param element The name of the element
    /// @param amount The same mole amount for the all elements (in units of mol)
    auto setElementAmount(std::string element, double amount) -> EquilibriumProblem&;

    /// Set the mole amount of electrical charge.
    /// @param amount The mole amount of electrical charge (in units of mol)
    auto setElectricalCharge(double amount) -> EquilibriumProblem&;

    /// Set the mole amount of an specie for the equilibrium calculation (in units of mol)
    /// @param species The reference for the specie in ChemicalSystem
    /// @param molar_amount The same mole amount for the all species (in units of mol)
    auto setSpecieAmount(const Species& species, double molar_amount) -> EquilibriumProblem&;

    /// Set the mole amount of an specie for the equilibrium calculation (in units of mol)
    /// @param ispecie The index of the specie
    /// @param molar_amount The same mole amount for the all species (in units of mol)
    auto setSpecieAmount(Index ispecie, double molar_amount) -> EquilibriumProblem&;

    /// Add a given amount of a compound or species to the equilibrium recipe.
    /// This method will first check if the given compound is present in the chemical system.
    /// If true, then `addSpecies` will be called. Otherwise, `addCompound` will be called.
    /// @param name The name of the compound or species
    /// @param amount The amount of the compound or species
    /// @param units The units of the amount (must be convertible to either mol or kg)
    auto add(std::string name, double amount, std::string units) -> EquilibriumProblem&;

    /// Add the mole amounts of the equilibrium species in a ChemicalState instance to the equilibrium recipe.
    /// This method only extract the mole amounts of equilibrium species in the given chemical state.
    /// If not all species in the system is under equilibrium assumption for this equilibrium calculation,
    /// ensure you set the partition of the chemical system before calling this method.
    /// @note If a multiplication factor is needed, for example `2.0`, use `add(2.0*state)`.
    /// @param state The ChemicalState instance with the mole amounts of the species
    /// @see setPartition
    auto add(const ChemicalState& state) -> EquilibriumProblem&;

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

    /// Add the mole amounts of the equilibrium species in a ChemicalState instance to the equilibrium recipe.
    /// This method only extract the mole amounts of equilibrium species in the given chemical state.
    /// If not all species in the system is under equilibrium assumption for this equilibrium calculation,
    /// ensure you set the partition of the chemical system before calling this method.
    /// @note If a multiplication factor is needed, for example `2.0`, use `addState(2.0*state)`.
    /// @param state The ChemicalState instance with the mole amounts of the species
    /// @see setPartition
    auto addState(const ChemicalState& state) -> EquilibriumProblem&;

    /// Return a reference to the ChemicalSystem instance used to create this EquilibriumProblem instance
    auto system() const -> const ChemicalSystem&;

    /// Return a reference to the Partition instance used to create this EquilibriumProblem instance
    auto partition() const -> const Partition&;

    /// Return the temperature for the equilibrium calculation (in units of K)
    auto temperature() const -> double;

    /// Return the pressure for the equilibrium calculation (in units of Pa)
    auto pressure() const -> double;

    /// Return the amounts of the elements for the equilibrium calculation (in units of mol)
    auto elementAmounts() const -> VectorConstRef;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
