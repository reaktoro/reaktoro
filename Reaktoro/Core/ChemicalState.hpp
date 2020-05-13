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

// C++ includes
#include <memory>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/PhaseChemicalProps.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalSystem;

/// The chemical state of a chemical system.
/// @see ChemicalSystem
/// @ingroup Core
class ChemicalState
{
public:
    /// Construct a ChemicalState instance with standard conditions.
    /// This constructor creates an instance of ChemicalState with temperature
    /// 25 °C, pressure 1 bar, and zero mole amounts for the species.
    explicit ChemicalState(const ChemicalSystem& system);

    /// Construct a copy of a ChemicalState instance.
    ChemicalState(const ChemicalState& other);

    /// Destroy this ChemicalState instance.
    virtual ~ChemicalState();

    /// Assign a ChemicalState instance to this instance.
    auto operator=(ChemicalState other) -> ChemicalState&;

    /// Set the temperature of the chemical state (in K).
    auto setTemperature(real val) -> void;

    /// Set the temperature of the chemical state with given unit.
    auto setTemperature(real val, String unit) -> void;

    /// Set the pressure of the chemical state (in Pa).
    auto setPressure(real val) -> void;

    /// Set the pressure of the chemical state with given unit.
    auto setPressure(real val, String unit) -> void;

    /// Set the amounts of all species in the chemical state to a common value (in mol).
    auto setSpeciesAmounts(real val) -> void;

    /// Set the amounts of the species in the chemical state (in mol).
    auto setSpeciesAmounts(ArrayXrConstRef n) -> void;

    /// Set the amount of a species with given index (in mol).
    auto setSpeciesAmount(Index ispecies, real amount) -> void;

    /// Set the amount of a species with given index and molar unit (convertible to mol).
    auto setSpeciesAmount(Index ispecies, real amount, String unit) -> void;

    /// Set the amount of a species with given name (in mol).
    auto setSpeciesAmount(String name, real amount) -> void;

    /// Set the amount of a species with given name and molar unit (convertible to mol).
    auto setSpeciesAmount(String name, real amount, String unit) -> void;

    /// Set the mass of a species with given index (in kg).
    auto setSpeciesMass(Index ispecies, real mass) -> void;

    /// Set the mass of a species with given index and mass unit (convertible to kg).
    auto setSpeciesMass(Index ispecies, real mass, String unit) -> void;

    /// Set the mass of a species with given name (in kg).
    auto setSpeciesMass(String name, real mass) -> void;

    /// Set the mass of a species with given name and mass unit (convertible to kg).
    auto setSpeciesMass(String name, real mass, String unit) -> void;

    /// Return the underlying chemical system for this chemical state.
    auto system() const -> const ChemicalSystem&;

    /// Return the temperature in the chemical state (in K).
    auto temperature() const -> real;

    /// Return the pressure in the chemical state (in Pa).
    auto pressure() const -> real;

    /// Return the amounts of the species in the chemical state (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the amounts of the elements in the chemical state (in mol).
    auto elementAmounts() const -> ArrayXr;

    /// Return the amount of the species in the chemical state with given index (in mol).
    auto speciesAmount(Index ispecies) const -> real;

    /// Return the amount of the species in the chemical state with given index and unit (convertible to mol).
    auto speciesAmount(Index ispecies, String unit) const -> real;

    /// Return the amount of the species in the chemical state with given name (in mol).
    auto speciesAmount(String name) const -> real;

    /// Return the amount of the species in the chemical state with given name and unit (convertible to mol).
    auto speciesAmount(String name, String unit) const -> real;

    /// Return the mass of the species in the chemical state with given index (in kg).
    auto speciesMass(Index ispecies) const -> real;

    /// Return the mass of the species in the chemical state with given index and unit (convertible to kg).
    auto speciesMass(Index ispecies, String unit) const -> real;

    /// Return the mass of the species in the chemical state with given name (in kg).
    auto speciesMass(String name) const -> real;

    /// Return the mass of the species in the chemical state with given name and unit (convertible to kg).
    auto speciesMass(String name, String unit) const -> real;

    /// Return the chemical properties of a phase in the chemical system.
    auto phaseProps(Index iphase) const -> PhaseChemicalProps;

    /// Return the chemical properties of the chemical system.
    auto props() const -> ChemicalProps;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
