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

namespace Reaktoro {

// Forward declarations
class PhaseChemicalProperties;
class ThermoProperties;

/// Defines the result of the function that calculates the chemical properties of a phase.
struct PhaseMixingModelResult
{
    /// The natural log of the activity constants of the species.
    ChemicalVector ln_activity_constants;

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The residual molar Gibbs energy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_gibbs_energy;

    /// The residual molar enthalpy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_enthalpy;

    /// The residual molar volume of the phase w.r.t. to its ideal state (in units of m3/mol).
    ChemicalScalar residual_molar_volume;

    /// The residual molar isobaric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cp;

    /// The residual molar isochoric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cv;
};

/// Defines the function signature for the calculation of chemical properties of a phase.
using PhaseMixingModel = std::function<PhaseMixingModelResult(double, double, const Vector&)>;

/// Defines the enumeration of possible standard reference states for a phase.
enum PhaseReferenceStateType { IdealGas, IdealSolution };

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

    /// Set the function that calculates the thermodynamic properties of the phase.
    auto setMixingModel(const PhaseMixingModel& model) -> void;

    /// Set the standard reference state type of the phase (default: IdealGas).
    auto setReferenceStateType(PhaseReferenceStateType reftype) -> void;

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

    /// Return the standard reference state type of the phase.
    auto referenceStateType() const -> PhaseReferenceStateType;

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
