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
#include <Reaktoro/Core/ChemicalSystem.hpp>

// Forward declarations (Optima)
namespace Optima { class State; }

namespace Reaktoro {

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
    explicit ChemicalState(ChemicalSystem const& system);

    /// Construct a copy of a ChemicalState instance.
    ChemicalState(ChemicalState const& other);

    /// Destroy this ChemicalState instance.
    virtual ~ChemicalState();

    /// Assign a ChemicalState instance to this instance.
    auto operator=(ChemicalState other) -> ChemicalState&;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING/GETTING TEMPERATURE
    // --------------------------------------------------------------------------------------------

    /// Set the temperature of the chemical state.
    /// @param value The temperature value (in K).
    auto setTemperature(real const& value) -> void;

    /// Set the temperature of the chemical state with given unit.
    /// @param value The temperature value.
    /// @param unit The temperature unit (must be convertible to K).
    auto setTemperature(real value, Chars unit) -> void;

    /// Set the temperature of the chemical state.
    /// @param value The temperature value (in K).
    /// @note This method is equivalent to ChemicalState::setTemperature(real)
    auto temperature(real const& value) -> void;

    /// Set the temperature of the chemical state with given unit.
    /// @param value The temperature value.
    /// @param unit The temperature unit (must be convertible to K).
    /// @note This method is equivalent to ChemicalState::setTemperature(real, String)
    auto temperature(real value, Chars unit) -> void;

    /// Return the temperature in the chemical state (in K).
    auto temperature() const -> real;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING/GETTING PRESSURE
    // --------------------------------------------------------------------------------------------

    /// Set the pressure of the chemical state.
    /// @param value The pressure value (in Pa).
    auto setPressure(real const& value) -> void;

    /// Set the pressure of the chemical state with given unit.
    /// @param value The pressure value.
    /// @param unit The pressure unit (must be convertible to Pa).
    auto setPressure(real value, Chars unit) -> void;

    /// Set the pressure of the chemical state.
    /// @param value The pressure value (in Pa).
    /// @note This method is equivalent to ChemicalState::setPressure(real)
    auto pressure(real const& value) -> void;

    /// Set the pressure of the chemical state with given unit.
    /// @param value The pressure value.
    /// @param unit The pressure unit (must be convertible to Pa).
    /// @note This method is equivalent to ChemicalState::setPressure(real, String)
    auto pressure(real value, Chars unit) -> void;

    /// Return the pressure in the chemical state (in Pa).
    auto pressure() const -> real;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING THE AMOUNT OR MASS OF SPECIES
    // --------------------------------------------------------------------------------------------

    /// Set the amounts of all species in the chemical state to a common value (in mol).
    auto setSpeciesAmounts(real const& value) -> void;

    /// Set the amounts of the species in the chemical state with given array (in mol).
    auto setSpeciesAmounts(ArrayXrConstRef const& n) -> void;

    /// Set the amounts of the species in the chemical state with given array (in mol).
    auto setSpeciesAmounts(ArrayXdConstRef const& n) -> void;

    /// Set the amount of a specific species in the system (in mol).
    /// @param ispecies The index of the species.
    /// @param amount The amount of the species
    auto setSpeciesAmount(Index ispecies, real const& amount) -> void;

    /// Set the amount of a specific species in the system.
    /// @param species The name or index of the species.
    /// @param amount The amount of the species
    /// @param unit The amount unit (must be convertible to mol)
    auto setSpeciesAmount(StringOrIndex const& species, real amount, Chars unit) -> void;

    /// Set the mass of a specific species in the system.
    /// @param species The name or index of the species.
    /// @param mass The mass of the species
    /// @param unit The mass unit (must be convertible to kg)
    auto setSpeciesMass(StringOrIndex const& species, real mass, Chars unit) -> void;

    /// Set the amount or mass of a chemical species in the chemical state.
    /// @param species The name or index of the species.
    /// @param value The amount or mass value of the species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if the chemical system has no species with name `species`.
    auto set(StringOrIndex const& species, real value, Chars unit) -> void;

    /// Add a specified amount or mass of a chemical species in the chemical state.
    /// @param species The name or index of the species in the chemical system.
    /// @param value The amount or mass value of the added species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if the chemical system has no species with name `species`.
    auto add(StringOrIndex const& species, real value, Chars unit) -> void;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR GETTING THE AMOUNT OR MASS OF SPECIES, ELEMENTS, AND CHARGE
    // --------------------------------------------------------------------------------------------

    /// Return the amounts of the species in the chemical state (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef;

    /// Return the amounts of the species from a phase in the chemical state (in mol).
    /// @param phase The name or index of the phase in the system.
    auto speciesAmountsInPhase(StringOrIndex const& phase) const -> ArrayXrConstRef;

    /// Return the amount of a species in the chemical state (in mol).
    /// @param species The name or index of the species.
    auto speciesAmount(StringOrIndex const& species) const -> real;

    /// Return the mass of a species in the chemical state (in kg).
    /// @param species The name or index of the species.
    auto speciesMass(StringOrIndex const& species) const -> real;

    /// Return the amounts of the conservative components (elements and charge) in the chemical state (in mol).
    auto componentAmounts() const -> ArrayXr;

    /// Return the amounts of the elements in the chemical state (in mol).
    auto elementAmounts() const -> ArrayXr;

    /// Return the electric charge in the chemical state (in mol).
    auto charge() const -> real;

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE AMOUNTS OF SPECIES IN THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    /// Scale the amounts of every species by a given factor.
    /// @param scalar The scale factor
    auto scaleSpeciesAmounts(real const& scalar) -> void;

    /// Scale the amounts of the species with given indices.
    /// @param scalar The scale factor
    /// @param indices The indices of the species
    auto scaleSpeciesAmounts(real const& scalar, Indices const& indices) -> void;

    /// Scale the amounts of the species in a phase by a given factor.
    /// @param phase The name or index of the phase in the system.
    /// @param scalar The scale factor
    auto scaleSpeciesAmountsInPhase(StringOrIndex const& phase, real const& scalar) -> void;

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE VOLUME OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    /// Scale the volume of the system to a new value.
    /// @param volume The new volume of the system
    /// @param unit The volume unit (must be convertible to m3)
    auto scaleVolume(real volume, Chars unit) -> void;

    /// Scale the volume of a phase to a new value.
    /// @param phase The name or index of the phase in the system.
    /// @param volume The new volume of the phase
    /// @param unit The volume unit (must be convertible to m3)
    auto scalePhaseVolume(StringOrIndex const& phase, real volume, Chars unit) -> void;

    /// Scale the total volume of fluids in the system to a new value.
    /// @param volume The new total volume of fluids in the system
    /// @param unit The volume unit (must be convertible to m3)
    auto scaleFluidVolume(real volume, Chars unit) -> void;

    /// Scale the total volume of solids in the system to a new value.
    /// @param volume The new total volume of solids in the system
    /// @param unit The volume unit (must be convertible to m3)
    auto scaleSolidVolume(real volume, Chars unit) -> void;

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE MASS OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    /// Scale the mass of the system to a new value.
    /// @param mass The new mass of the system
    /// @param unit The mass unit (must be convertible to kg)
    auto scaleMass(real mass, Chars unit) -> void;

    /// Scale the mass of a phase to a new value.
    /// @param phase The name or index of the phase in the system.
    /// @param mass The new mass of the phase
    /// @param unit The mass unit (must be convertible to kg)
    auto scalePhaseMass(StringOrIndex const& phase, real mass, Chars unit) -> void;

    /// Scale the total mass of fluids in the system to a new value.
    /// @param mass The new total mass of fluids in the system
    /// @param unit The mass unit (must be convertible to kg)
    auto scaleFluidMass(real mass, Chars unit) -> void;

    /// Scale the total mass of solids in the system to a new value.
    /// @param mass The new total mass of solids in the system
    /// @param unit The mass unit (must be convertible to kg)
    auto scaleSolidMass(real mass, Chars unit) -> void;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING/GETTING SURFACE AREAS BETWEEN PHASES
    // --------------------------------------------------------------------------------------------

    /// Set the surface area between two phases.
    /// @param s The surface areas of reacting phase interfaces in the system (in m2)
    auto setSurfaceAreas(ArrayXrConstRef const& s) -> void;

    /// Set the surface area between two phases.
    /// @param s The surface areas of reacting phase interfaces in the system (in m2)
    auto setSurfaceAreas(ArrayXdConstRef const& s) -> void;

    /// Set the surface area between two phases.
    /// @param phase1 The name or index of a phase.
    /// @param phase2 The name or index of the other phase.
    /// @param value The surface area between these two phases.
    /// @param unit The unit of the surface area (must be convertible to m2).
    auto setSurfaceArea(StringOrIndex const& phase1, StringOrIndex const& phase2, real value, Chars unit) -> void;

    /// Set the surface area between two phases with given surface index.
    /// @param isurface The index of the surface.
    /// @param value The surface area between these two phases.
    /// @param unit The unit of the surface area (must be convertible to m2).
    auto setSurfaceArea(Index isurface, real value, Chars unit) -> void;

    /// Set the surface area between two phases.
    /// @param phase1 The name or index of a phase.
    /// @param phase2 The name or index of the phase interfacing with the previous one.
    /// @param value The surface area value.
    /// @param unit The unit of the surface area (must be convertible to m2).
    /// @note This method is equivalent to ChemicalState::setSurfaceArea
    auto surfaceArea(StringOrIndex const& phase1, StringOrIndex const& phase2, real value, Chars unit) -> void;

    /// Return the surface area of the interface between two reacting phases (in m2).
    /// @param phase1 The name or index of a phase.
    /// @param phase2 The name or index of the phase interfacing with the previous one.
    /// @warning An error is thrown if no surface area has been set for the phase pair `phase1` and `phase2`.
    auto surfaceArea(StringOrIndex const& phase1, StringOrIndex const& phase2) const -> real;

    /// Return the surface area of the interface between two reacting phases with given surface index (in m2).
    /// @param isurface The index of the surface between two phases.
    auto surfaceArea(Index isurface) const -> real;

    /// Return the areas of all reacting phase interfaces in the system (in m2).
    auto surfaceAreas() const -> ArrayXrConstRef;

    /// Return the phase pairs for which the surface area has been defined.
    /// @see surfaceArea, setSurfaceArea
    auto surfaces() const -> const Pairs<Index, Index>&;

    // --------------------------------------------------------------------------------------------
    // METHODS FOR UPDATING CHEMICAL STATE AND ITS PROPERTIES
    // --------------------------------------------------------------------------------------------

    /// Update the chemical state and properties of the system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    auto update(real const& T, real const& P, ArrayXrConstRef const& n) -> void;

    /// Update the chemical state and properties of the system.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    /// @param s The surface areas of reacting phase interfaces in the system (in m2)
    auto update(real const& T, real const& P, ArrayXrConstRef const& n, ArrayXrConstRef const& s) -> void;

    /// Update the chemical state and properties of the system using ideal activity models.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    auto updateIdeal(real const& T, real const& P, ArrayXrConstRef const& n) -> void;

    /// Update the chemical state and properties of the system using ideal activity models.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    /// @param s The surface areas of reacting phase interfaces in the system (in m2)
    auto updateIdeal(real const& T, real const& P, ArrayXrConstRef const& n, ArrayXrConstRef const& s) -> void;

    // --------------------------------------------------------------------------------------------
    // MISCELLANEOUS METHODS
    // --------------------------------------------------------------------------------------------

    /// Return the underlying chemical system for this chemical state.
    auto system() const -> ChemicalSystem const&;

    /// Return the chemical properties of the system. For performance reasons,
    /// the stored chemical properties are not updated at every change in the
    /// chemical state. For a ChemicalState object `state`, update its chemical
    /// properties using `state.props().update(state)`.
    auto props() const -> ChemicalProps const&;

    /// Return the chemical properties of the system. For performance reasons,
    /// the stored chemical properties are not updated at every change in the
    /// chemical state. For a ChemicalState object `state`, update its chemical
    /// properties using `state.props().update(state)`.
    auto props() -> ChemicalProps&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() const -> Equilibrium const&;

    /// Return the equilibrium properties of a calculated chemical equilibrium state.
    auto equilibrium() -> Equilibrium&;

    /// Output this ChemicalState instance to a stream.
    auto output(std::ostream& out) const -> void;

    /// Output this ChemicalState instance to a file.
    auto output(String const& filename) const -> void;

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
    Equilibrium(ChemicalSystem const& system);

    /// Construct a copy of a ChemicalState::Equilibrium instance
    Equilibrium(Equilibrium const& other);

    /// Destroy this ChemicalState::Equilibrium instance
    virtual ~Equilibrium();

    /// Assign a ChemicalState::Equilibrium instance to this instance
    auto operator=(Equilibrium other) -> Equilibrium&;

    /// Set the names of the input variables used in the equilibrium calculation.
    auto setInputNames(Strings const& names) -> void;

    /// Set the values of the input variables used in the equilibrium calculation.
    auto setInputValues(ArrayXdConstRef const& w) -> void;

    /// Set initial component amounts used in the equilibrium calculation.
    auto setInitialComponentAmounts(ArrayXdConstRef const& c0) -> void;

    /// Set the computed control variables *p* in the equilibrium calculation.
    auto setControlVariablesP(ArrayXdConstRef const& p) -> void;

    /// Set the computed control variables *q* in the equilibrium calculation.
    auto setControlVariablesQ(ArrayXdConstRef const& q) -> void;

    /// Set the Optima::State object computed as part of the equilibrium calculation.
    auto setOptimaState(Optima::State const& state) -> void;

    /// Return the number of primary species.
    auto numPrimarySpecies() const -> Index;

    /// Return the number of secondary species.
    auto numSecondarySpecies() const -> Index;

    /// Return the indices of the primary species.
    auto indicesPrimarySpecies() const -> ArrayXlConstRef;

    /// Return the indices of the secondary species.
    auto indicesSecondarySpecies() const -> ArrayXlConstRef;

    /// Return the chemical potentials of the elements in the equilibrium state (in unit of J/mol).
    auto elementChemicalPotentials() const -> ArrayXdConstRef;

    /// Return the stabilities of the chemical species in the equilibrium state (normalized by *RT*).
    /// These are the slack variables with respect to bound constraints on the
    /// amounts of the species in a chemical equilibrium calculation. They can
    /// be interpreted as measures of stability of a species at equilibrium,
    /// with values closer to zero meaning more stable and away from the bounds.
    auto speciesStabilities() const -> ArrayXdConstRef;

    /// Return the amounts of the explicit titrants in the equilibrium state (in mol).
    auto explicitTitrantAmounts() const -> ArrayXdConstRef;

    /// Return the amounts of the implicit titrants in the equilibrium state (in mol).
    auto implicitTitrantAmounts() const -> ArrayXdConstRef;

    /// Return the names of the input variables used in the equilibrium calculation.
    auto inputNames() const -> Strings const&;

    /// Return the values of the input variables used in the equilibrium calculation.
    auto inputValues() const -> ArrayXdConstRef;

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
    auto w() const -> ArrayXdConstRef;

    /// Return the initial component amounts used in the equilibrium calculation.
    auto c() const -> ArrayXdConstRef;

    /// Return the Optima::State object computed as part of the equilibrium calculation.
    auto optimaState() const -> Optima::State const&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

/// Output a ChemicalState object to an output stream.
auto operator<<(std::ostream& out, ChemicalState const& state) -> std::ostream&;

} // namespace Reaktoro

//=========================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING CHEMICALSTATE
//=========================================================================

namespace Reaktoro {

template<typename T>
struct MemoizationTraits;

/// Specialize MemoizationTraits for ChemicalState.
template<>
struct MemoizationTraits<ChemicalState>
{
    using Type = ChemicalState;

    /// The type used instead to cache a ChemicalState object.
    using CacheType = Tuple<real, real, ArrayXr, ArrayXr>;

    auto equal(const Tuple<real, real, ArrayXr, ArrayXr>& a, ChemicalState const& b)
    {
        auto const& [T, P, n, S] = a;
        return
            T == b.temperature() &&
            P == b.pressure() &&
            (n == b.speciesAmounts()).all() &&
            (S == b.surfaceAreas()).all();
    }

    auto assign(Tuple<real, real, ArrayXr, ArrayXr>& a, ChemicalState const& b)
    {
        auto& [T, P, n, S] = a;
        T = b.temperature();
        P = b.pressure();
        n = b.speciesAmounts();
        S = b.surfaceAreas();
    }
};

} // namespace Reaktoro
