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

/// A type used to define a phase and its attributes.
/// @see ChemicalSystem, Element, Species
/// @ingroup Core
class Phase
{
public:
    /// Construct a default Phase instance.
    Phase();

    /// Construct a copy of a Phase instance
    Phase(const Phase& other);

    /// Destroy this instance
    virtual ~Phase();

    /// Assign an Phase instance to this instance
    auto operator=(Phase other) -> Phase&;

    /// Set the name of the phase.
    auto setName(std::string name) -> void;

    /// Set the species of the phase.
    auto setSpecies(const std::vector<Species>& species) -> void;

    /// Set the function for the concentrations of the species in the phase.
    auto setConcentrationFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the natural log of the activity coefficients of the species in the phase.
    auto setActivityCoefficientFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the natural log of the activities of the species in the phase.
    auto setActivityFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the molar volume of the phase (in units of m3/mol).
    /// If the molar volume function of the phase is not set, then a default function
    /// based on the standard molar volumes of the species will be used:
    /// @f[v_{\pi}=\sum_{i}x_{i}v_{i}^{\circ},@f]
    /// where @f$v_{\pi}@f$ is the molar volume of the phase; @f$x_{i}@f$ and
    /// @f$v_{i}@f$ are the molar fraction and standard molar volume of the *i*-th species.
    auto setMolarVolumeFunction(const ChemicalScalarFunction& function) -> void;

    /// Return the number of elements in the phase
    auto numElements() const -> unsigned;

    /// Return the number of species in the phase
    auto numSpecies() const -> unsigned;

    /// Return the name of the phase
    auto name() const -> std::string;

    /// Return the elements of the phase
    auto elements() const -> const std::vector<Element>&;

    /// Return the species of the phase
    auto species() const -> const std::vector<Species>&;

    /// Return the species of the phase with a given index
    /// @param index The index of the species
    auto species(Index index) const -> const Species&;

    /// Return the function for the concentrations of the species in the phase.
    auto concentrationFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the natural log of the activity coefficients of the species in the phase.
    auto activityCoefficientFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the natural log of the activities of the species in the phase.
    auto activityFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the molar volume of the phase (in units of m3/mol).
    auto molarVolumeFunction() const -> const ChemicalScalarFunction&;

    /// Calculate the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    auto standardGibbsEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar enthalpies of the species (in units of J/mol).
    auto standardEnthalpies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    auto standardHelmholtzEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar entropies of the species (in units of J/K).
    auto standardEntropies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar volumes of the species (in units of m3/mol).
    auto standardVolumes(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar internal energies of the species (in units of J/mol).
    auto standardInternalEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto standardHeatCapacities(double T, double P) const -> ThermoVector;

    /// Calculate the concentrations of the species (no uniform units).
    auto concentrations(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activity coefficients of the species.
    auto activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activities of the species.
    auto activities(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the molar volume of the phase (in units of m3/mol).
    auto molarVolume(double T, double P, const Vector& n) const -> ChemicalScalar;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

/// Compare two Phase instances for less than
auto operator<(const Phase& lhs, const Phase& rhs) -> bool;

/// Compare two Phase instances for equality
auto operator==(const Phase& lhs, const Phase& rhs) -> bool;

} // namespace Reaktoro
