// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

// Forward declarations (Optima)
namespace Optima { class State; }

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class Params;

//=================================================================================================
//
// ChemicalState
//
//=================================================================================================

/// The chemical state of a chemical system.
/// @see ChemicalSystem
/// @ingroup Core
class ChemicalState
{
public:
    // Forward declarations
    class Equilibrium;
    class Props;

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

    /// Set the amounts of the species in the chemical state (in mol).
    auto setSpeciesAmounts(ArrayXdConstRef n) -> void;

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

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() const -> const Equilibrium&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() -> Equilibrium&;

    /// Return the chemical properties of the system.
    /// @warning For performance reasons, the stored chemical properties are
    /// not updated at every change in the chemical state. For a ChemicalState
    /// object `state`, update its chemical properties using `state.props().update()`.
    auto props() const -> const Props&;

    /// Return the chemical properties of the system.
    /// @warning For performance reasons, the stored chemical properties are
    /// not updated at every change in the chemical state. For a ChemicalState
    /// object `state`, update its chemical properties using `state.props().update()`.
    auto props() -> Props&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

//=================================================================================================
//
// ChemicalState::Equilibrium
//
//=================================================================================================

/// The access to the properties related to an equilibrium state in a ChemicalState object.
class ChemicalState::Equilibrium
{
public:
    /// Construct a ChemicalState::Equilibrium instance.
    Equilibrium(const ChemicalSystem& system);

    /// Construct a copy of a ChemicalState::Equilibrium instance
    Equilibrium(const Equilibrium& other);

    /// Destroy this ChemicalState::Equilibrium instance
    virtual ~Equilibrium();

    /// Assign a ChemicalState::Equilibrium instance to this instance
    auto operator=(Equilibrium other) -> Equilibrium&;

    /// Set the input parameters representing the conditions used for the calculation of the equilibrium state.
    auto setParams(const Params& w) -> void;

    /// Set initial component amounts for the calculation of the equilibrium state.
    auto setInitialComponentAmounts(ArrayXdConstRef b) -> void;

    /// Set the Optima::State object obtained during an Optima optimization calculation.
    auto setOptimaState(const Optima::State& state) -> void;

    /// Return the last set Optima::State object.
    auto optimaState() const -> const Optima::State&;

    /// Set the indices of the primary and secondary species at the equilibrium state.
    /// @param ips The indices of the equilibrium species ordered as (primary, secondary)
    /// @param kp The number of primary species
    auto setIndicesPrimarySecondarySpecies(ArrayXlConstRef ips, Index kp) -> void;

    /// Set the indices of elements whose amounts should be positive, but given amount was less or equal to zero.
    auto setIndicesStrictlyUnstableElements(ArrayXlConstRef isue) -> void;

    /// Set the indices of species that contain one or more strictly unstable elements.
    /// @see setIndicesElementsStrictlyUnstable
    auto setIndicesStrictlyUnstableSpecies(ArrayXlConstRef isus) -> void;

    /// Return the number of primary species.
    auto numPrimarySpecies() const -> Index;

    /// Return the number of secondary species.
    auto numSecondarySpecies() const -> Index;

    /// Return the indices of the primary species.
    auto indicesPrimarySpecies() const -> ArrayXlConstRef;

    /// Return the indices of the secondary species.
    auto indicesSecondarySpecies() const -> ArrayXlConstRef;

    /// Return the indices of elements whose amounts should be positive, but given amount was less or equal to zero.
    auto indicesStrictlyUnstableElements() const -> ArrayXlConstRef;

    /// Return the indices of species that contain one or more strictly unstable elements.
    auto indicesStrictlyUnstableSpecies() const -> ArrayXlConstRef;

    /// Return the chemical potentials of the elements in the equilibrium state (in unit of J/mol).
    auto elementChemicalPotentials() const -> ArrayXdConstRef;

    /// Return the stabilities of the chemical species in the equilibrium state (in unit of J/mol).
    /// These are the slack variables with respect to bound constraints on the
    /// amounts of the species in a chemical equilibrium calculation. They can
    /// be interpreted as measures of stability of a species at equilibrium,
    /// with values closer to zero meaning more stable and away from the bounds.
    auto speciesStabilities() const -> ArrayXdConstRef;

    /// Return the amounts of the explicit titrants in the equilibrium state (in unit of mol).
    auto explicitTitrantAmounts() const -> ArrayXdConstRef;

    /// Return the amounts of the implicit titrants in the equilibrium state (in unit of mol).
    auto implicitTitrantAmounts() const -> ArrayXdConstRef;

    /// Return the control variables *p* in the equilibrium state.
    auto p() const -> ArrayXdConstRef;

    /// Return the control variables *q* in the equilibrium state.
    auto q() const -> ArrayXdConstRef;

    /// Return the input parameters for the calculation of the equilibrium state.
    auto w() const -> const Params&;

    /// Return the initial component amounts for the calculation of the equilibrium state.
    auto b() const -> ArrayXdConstRef;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

//=================================================================================================
//
// ChemicalState::Props
//
//=================================================================================================

/// The access to the chemical properties of the system.
class ChemicalState::Props : public ChemicalProps
{
public:
    /// Construct a ChemicalState::Props instance with given associated chemical state.
    Props(const ChemicalSystem& system, const ChemicalState& state);

    /// Construct a copy of a ChemicalState::Props instance
    Props(const Props& other);

    /// Destroy this ChemicalState::Props instance
    virtual ~Props();

    /// Assign a ChemicalState::Props instance to this instance
    auto operator=(Props other) -> Props&;

    /// Update the chemical properties in the system according to its current state.
    auto update() -> void;

    // Inherit all other update methods.
    using ChemicalProps::update;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
