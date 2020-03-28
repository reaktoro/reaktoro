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
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The possible physical states of a phase.
/// @see Phase
enum class PhasePhysicalState
{
    Solid, Liquid, Gas, Plasma
};

/// A type used to define a phase and its attributes.
/// @see ChemicalSystem, Element, Species
/// @ingroup Core
class Phase
{
public:
    /// Construct a default Phase object.
    Phase();

    /// Construct a copy of a Phase object.
    Phase(const Phase& other);

    /// Destroy this Phase object.
    ~Phase();

    /// Assign another Phase object to this.
    auto operator=(Phase other) -> Phase&;

    /// Set the name of the phase.
    auto setName(std::string name) -> void;

    /// Set the type of the phase.
    auto setType(std::string type) -> void;

    /// Set the physical state of the phase.
    auto setPhysicalState(PhasePhysicalState state) -> void;

    /// Set the species of the phase.
    auto setSpecies(const std::vector<Species>& species) -> void;

    /// Set the standard thermodynamic model of the phase.
    auto setStandardThermoModel(const StandardThermoModelFn& model) -> void;

    /// Set the activity model of the phase.
    auto setActivityModel(const ActivityModelFn& model) -> void;

    /// Return the number of elements in the phase.
    auto numElements() const -> unsigned;

    /// Return the number of species in the phase.
    auto numSpecies() const -> unsigned;

    /// Return the name of the phase.
    auto name() const -> std::string;

    /// Return the type of the phase.
    auto type() const -> std::string;

    /// Return the physical state of the phase.
    auto physicalState() const -> PhasePhysicalState;

    /// Return the elements of the phase.
    auto elements() const -> const std::vector<Element>&;

    /// Return the species of the phase.
    auto species() const -> const std::vector<Species>&;

    /// Return the species of the phase with a given index.
    auto species(Index index) const -> const Species&;

    /// Return true if the state of matter of the phase is liquid, gas or plasma.
    auto isFluid() const -> bool;

    /// Return true if the phase type is solid.
    auto isSolid() const -> bool;

    /// Return the standard thermodynamic model of the phase.
    /// @see StandardThermoModel
    auto standardThermoModel() const -> const StandardThermoModelFn&;

    /// Return the activity model of the phase.
    /// @see ActivityModel
    auto activityModel() const -> const ActivityModelFn&;

    /// Return the index of a species in the phase.
    /// @param name The name of the species
    /// @return The index of the species if found, or the number of species in the phase otherwise.
    auto indexSpecies(std::string name) const -> Index;

    /// Return the index of a species in the system.
    /// @param name The name of the species
    /// @return The index of the species if found, or a runtime exception otherwise.
    auto indexSpeciesWithError(std::string name) const -> Index;

    /// Return the index of the first species in the phase with any of the given names.
    /// @param names The tentative names of the species in the phase.
    /// @return The index of the species if found, or the number of species in the phase otherwise.
    auto indexSpeciesAny(const std::vector<std::string>& names) const -> Index;

    /// Return the index of the first species in the phase with any of the given names.
    /// @param names The tentative names of the species in the phase.
    /// @return The index of the species if found, or a runtime exception otherwise.
    auto indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index;

    /// Calculate the standard thermodynamic properties of the phase.
    /// @param T The temperature of the system (in units of K)
    /// @param P The pressure of the system (in units of Pa)
    auto standardThermoProps(double T, double P) const -> StandardThermoProps;

    /// Calculate the activity and excess thermodynamic properties of the phase.
    /// @param T The temperature of the system (in units of K)
    /// @param P The pressure of the system (in units of Pa)
    /// @param n The amounts of the species in the phase (in units of mol)
    auto activityProps(double T, double P, VectorXrConstRef n) const -> ActivityProps;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Compare two Phase instances for less than
auto operator<(const Phase& lhs, const Phase& rhs) -> bool;

/// Compare two Phase instances for equality
auto operator==(const Phase& lhs, const Phase& rhs) -> bool;

} // namespace Reaktoro
