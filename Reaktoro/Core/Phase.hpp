// Reaktoro is a C++ library for computational reaction modelling.
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
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/PhaseThermoModel.hpp>

namespace Reaktoro {

// Forward declarations
class PhaseChemicalProperties;
class ThermoProperties;

/// Defines the enumeration of possible standard reference states for a phase.
enum class PhaseReferenceState { IdealGas, IdealSolution };

/// A type used to define a phase and its attributes.
/// @see ChemicalSystem, Element, Species
/// @ingroup Core
class Phase
{
public:
    /// Construct a default Phase instance.
    Phase();

    /// Construct a copy of a Phase instance.
    Phase(const Phase& other);

    /// Destroy this instance.
    virtual ~Phase();

    /// Assign an Phase instance to this instance.
    auto operator=(Phase other) -> Phase&;

    /// Set the name of the phase.
    auto setName(std::string name) -> void;

    /// Set the species of the phase.
    auto setSpecies(const std::vector<Species>& species) -> void;

    /// Set the standard reference state type of the phase (default: IdealGas).
    auto setReferenceState(PhaseReferenceState reftype) -> void;

    /// Set the function that calculates the standard thermodynamic properties of the phase.
    auto setThermoModel(const PhaseThermoModel& model) -> void;

    /// Set the function that calculates the chemical properties of the phase.
    auto setChemicalModel(const PhaseChemicalModel& model) -> void;

    /// Return the number of elements in the phase.
    auto numElements() const -> unsigned;

    /// Return the number of species in the phase.
    auto numSpecies() const -> unsigned;

    /// Return the name of the phase.
    auto name() const -> std::string;

    /// Return the elements of the phase.
    auto elements() const -> const std::vector<Element>&;

    /// Return the species of the phase.
    auto species() const -> const std::vector<Species>&;

    /// Return the species of the phase with a given index.
    auto species(Index index) const -> const Species&;

    /// Return the thermodynamic model function of the phase.
    auto thermoModel() const -> const PhaseThermoModel&;

    /// Return the thermodynamic model function of the phase.
    auto chemicalModel() const -> const PhaseChemicalModel&;

    /// Return the index of a species in the phase.
    auto indexSpecies(std::string species) const -> Index;

    /// Return the standard reference state type of the phase.
    auto referenceState() const -> PhaseReferenceState;

    /// Return the calculated standard thermodynamic properties of the species.
    auto properties(double T, double P) const -> ThermoProperties;

    /// Return the calculated chemical properties of the phase and its species.
    auto properties(double T, double P, const Vector& n) const -> PhaseChemicalProperties;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Compare two Phase instances for less than
auto operator<(const Phase& lhs, const Phase& rhs) -> bool;

/// Compare two Phase instances for equality
auto operator==(const Phase& lhs, const Phase& rhs) -> bool;

} // namespace Reaktoro
