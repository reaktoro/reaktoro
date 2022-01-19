// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

    /// Set the temperature of the chemical state.
    /// @param value The temperature value (in K).
    auto temperature(real value) -> void;

    /// Set the temperature of the chemical state with given unit.
    /// @param value The temperature value.
    /// @param unit The temperature unit (convertible to K).
    auto temperature(real value, String unit) -> void;

    /// Set the pressure of the chemical state.
    /// @param value The pressure value (in Pa).
    auto pressure(real value) -> void;

    /// Set the pressure of the chemical state with given unit.
    /// @param value The pressure value.
    /// @param unit The pressure unit (convertible to Pa).
    auto pressure(real value, String unit) -> void;

    /// Add a specified amount or mass of a chemical species in the chemical state.
    /// @param species The name of the species in the chemical system.
    /// @param value The amount or mass value of the added species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if the chemical system has no species with name @p species.
    auto add(String species, real value, String unit) -> void;

    /// Add a specified amount or mass of a chemical species in the chemical state.
    /// @param ispecies The index of the species in the chemical system.
    /// @param value The amount or mass value of the added species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto add(Index ispecies, real value, String unit) -> void;

    /// Set a specified amount or mass of a chemical species in the chemical state.
    /// @param species The name of the species in the chemical system.
    /// @param value The amount or mass value of the added species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if the chemical system has no species with name @p species.
    auto set(String species, real value, String unit) -> void;

    /// Set a specified amount or mass of a chemical species in the chemical state.
    /// @param ispecies The index of the species in the chemical system.
    /// @param value The amount or mass value of the added species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto set(Index ispecies, real value, String unit) -> void;

    /// Set the temperature of the chemical state.
    /// @param value The temperature value (in K).
    /// @note This method is equivalent to ChemicalState::temperature(real)
    auto setTemperature(real value) -> void;

    /// Set the temperature of the chemical state with given unit.
    /// @param value The temperature value.
    /// @param unit The temperature unit (convertible to K).
    /// @note This method is equivalent to ChemicalState::temperature(real, String)
    auto setTemperature(real value, String unit) -> void;

    /// Set the pressure of the chemical state.
    /// @param value The pressure value (in Pa).
    /// @note This method is equivalent to ChemicalState::pressure(real)
    auto setPressure(real value) -> void;

    /// Set the pressure of the chemical state with given unit.
    /// @param value The pressure value.
    /// @param unit The pressure unit (convertible to Pa).
    /// @note This method is equivalent to ChemicalState::pressure(real, String)
    auto setPressure(real value, String unit) -> void;

    /// Set the amounts of all species in the chemical state to a common value (in mol).
    auto setSpeciesAmounts(real value) -> void;

    /// Set the amounts of the species in the chemical state (in mol).
    auto setSpeciesAmounts(ArrayXrConstRef n) -> void;

    /// Set the amounts of the species in the chemical state (in mol).
    auto setSpeciesAmounts(ArrayXdConstRef n) -> void;

    /// Set the amount of a species in the system.
    /// @param species The name or index of the species in the system.
    /// @param amount The new amount of the species
    /// @param unit The amount unit
    auto setSpeciesAmount(StringOrIndex species, real amount, String unit) -> void;

    /// Set the mass of a species in the system.
    /// @param species The name or index of the species in the system.
    /// @param mass The new mass of the species
    /// @param unit The mass unit
    auto setSpeciesMass(StringOrIndex species, real mass, String unit) -> void;

    /// Scale the amounts of every species by a given scalar.
    /// @param scalar The scale factor
    auto scaleSpeciesAmounts(real scalar) -> void;

    /// Scale the amounts of the species in a phase by a given scalar.
    /// @param phase The name or index of the phase in the system.
    /// @param scalar The scale factor
    auto scaleSpeciesAmountsInPhase(StringOrIndex phase, real scalar) -> void;

    /// Scale the volume of the system by adjusting the amounts of its species equally.
    /// @param value The new volume of the system
    /// @param unit The volume unit
    auto scaleVolume(real value, String unit) -> void;

    /// Scale the volume of a phase by adjusting the amounts of its species equally.
    /// @param phase The name or index of the phase in the system.
    /// @param value The new volume of the phase
    /// @param unit The volume unit
    auto scalePhaseVolume(StringOrIndex phase, real value, String unit) -> void;

    // /// Scale the total volume of fluids in the system by adjusting the volumes of fluid phases equally.
    // /// @param value The new total volume of fluids in the system
    // /// @param unit The volume unit
    // auto scaleFluidVolume(real value, String unit) -> void;

    // /// Scale the total volume of solids in the system by adjusting the volumes of solid phases equally.
    // /// @param value The new total volume of solids in the system
    // /// @param unit The volume unit
    // auto scaleSolidVolume(real value, String unit) -> void;

    /// Scale the mass of the system by adjusting the amounts of its species equally.
    /// @param value The new mass of the system
    /// @param unit The mass unit
    auto scaleMass(real value, String unit) -> void;

    /// Scale the mass of a phase by adjusting the amounts of its species equally.
    /// @param phase The name or index of the phase in the system.
    /// @param value The new mass of the phase
    /// @param unit The mass unit
    auto scalePhaseMass(StringOrIndex phase, real value, String unit) -> void;

    // /// Scale the total mass of fluids in the system by adjusting the masses of fluid phases equally.
    // /// @param value The new total mass of fluids in the system
    // /// @param unit The mass unit
    // auto scaleFluidMass(real value, String unit) -> void;

    // /// Scale the total mass of solids in the system by adjusting the masses of solid phases equally.
    // /// @param value The new total mass of solids in the system
    // /// @param unit The mass unit
    // auto scaleSolidMass(real value, String unit) -> void;

    /// Return the underlying chemical system for this chemical state.
    auto system() const -> const ChemicalSystem&;

    /// Return the temperature in the chemical state (in K).
    auto temperature() const -> real;

    /// Return the pressure in the chemical state (in Pa).
    auto pressure() const -> real;

    /// Return the amounts of the species in the chemical state (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the amounts of the species from a phase in the chemical state (in mol).
    /// @param phase The name or index of the phase in the system.
    auto speciesAmountsInPhase(StringOrIndex phase) const -> ArrayXrConstRef;

    /// Return the amounts of the elements and charge in the chemical state (in mol).
    auto componentAmounts() const -> ArrayXr;

    /// Return the amounts of the elements in the chemical state (in mol).
    auto elementAmounts() const -> ArrayXr;

    /// Return the electric charge in the chemical state (in mol).
    auto charge() const -> real;

    /// Return the amount of a species in the chemical state (in mol).
    /// @param species The name or index of the species in the system.
    auto speciesAmount(StringOrIndex species) const -> real;

    /// Return the mass of a species in the chemical state (in kg).
    /// @param species The name or index of the species in the system.
    auto speciesMass(StringOrIndex species) const -> real;

    /// Return the chemical properties of the system. For performance reasons,
    /// the stored chemical properties are not updated at every change in the
    /// chemical state. For a ChemicalState object `state`, update its chemical
    /// properties using `state.props().update(state)`.
    auto props() const -> const ChemicalProps&;

    /// Return the chemical properties of the system. For performance reasons,
    /// the stored chemical properties are not updated at every change in the
    /// chemical state. For a ChemicalState object `state`, update its chemical
    /// properties using `state.props().update(state)`.
    auto props() -> ChemicalProps&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() const -> const Equilibrium&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() -> Equilibrium&;

    /// Output this ChemicalState instance to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output this ChemicalState instance to a file.
    auto output(const String& filename) const -> void;

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

    /// Set the names of the input variables used in the equilibrium calculation.
    auto setInputNames(const Strings& names) -> void;

    /// Set the values of the input variables used in the equilibrium calculation.
    auto setInputValues(VectorXdConstRef w) -> void;

    /// Set initial component amounts used in the equilibrium calculation.
    auto setInitialComponentAmounts(ArrayXdConstRef b) -> void;

    /// Set the computed control variables *p* in the equilibrium calculation.
    auto setControlVariablesP(ArrayXdConstRef p) -> void;

    /// Set the computed control variables *q* in the equilibrium calculation.
    auto setControlVariablesQ(ArrayXdConstRef q) -> void;

    /// Set the Optima::State object computed as part of the equilibrium calculation.
    auto setOptimaState(const Optima::State& state) -> void;

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

    /// Return the names of the input variables used in the equilibrium calculation.
    auto inputNames() const -> const Strings&;

    /// Return the values of the input variables used in the equilibrium calculation.
    auto inputValues() const -> VectorXdConstRef;

    /// Return the initial component amounts used in the equilibrium calculation.
    auto initialComponentAmounts() const -> ArrayXdConstRef;

    /// Return the computed control variables *p* in the equilibrium calculation.
    auto controlVariablesP() const -> ArrayXdConstRef;

    /// Return the computed control variables *q* in the equilibrium calculation.
    auto controlVariablesQ() const -> ArrayXdConstRef;

    /// Return the control variables *p* computed in the equilibrium calculation.
    auto p() const -> ArrayXdConstRef;

    /// Return the control variables *q* computed in the equilibrium calculation.
    auto q() const -> ArrayXdConstRef;

    /// Return the values of the input variables used in the equilibrium calculation.
    auto w() const -> VectorXdConstRef;

    /// Return the initial component amounts used in the equilibrium calculation.
    auto b() const -> ArrayXdConstRef;

    /// Return the Optima::State object computed as part of the equilibrium calculation.
    auto optimaState() const -> const Optima::State&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output a ChemicalState object to an output stream.
auto operator<<(std::ostream& out, const ChemicalState& state) -> std::ostream&;

} // namespace Reaktoro
