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

// Reaktor includes
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Phase.hpp>

namespace Reaktor {

/// The type used to define the attributes and model functions of a ChemicalSystem instance
/// @see ChemicalSystem
/// @ingroup Core
struct ChemicalSystemData
{
    /// The list of phases in the chemical system
    PhaseList phases;

    /// The function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    ThermoVectorFunction gibbs_energies;

    /// The function for the apparent standard molar enthalpies of the species (in units of J/mol).
    ThermoVectorFunction enthalpies;

    /// The function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    ThermoVectorFunction helmholtz_energies;

    /// The function for the standard molar entropies of the species (in units of J/K).
    ThermoVectorFunction entropies;

    /// The function for the standard molar volumes of the species (in units of m3/mol).
    ThermoVectorFunction volumes;

    /// The function for the apparent standard molar internal energies of the species (in units of J/mol).
    ThermoVectorFunction internal_energies;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoVectorFunction heat_capacities_cp;

    /// The function for the concentrations of the species (no uniform units).
    ChemicalVectorFunction concentrations;

    /// The function for the natural log of the activity coefficients of the species.
    ChemicalVectorFunction ln_activity_coefficients;

    /// The function for the natural log of the activities of the species.
    ChemicalVectorFunction ln_activities;

    /// The function for the chemical potentials of the species (in units of J/mol).
    ChemicalVectorFunction chemical_potentials;

    /// The function for the densities of the phases (in units of kg/m3).
    ChemicalVectorFunction densities;
};

/// The type used to define a chemical system and its attributes
/// @see Species, Phase
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a default ChemicalSystem instance
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with all its attributes
    ChemicalSystem(const ChemicalSystemData& data);

    /// Get the list of elements in the chemical system
    auto elements() const -> const ElementList&;

    /// Get the list of species in the chemical system
    auto species() const -> const SpeciesList&;

    /// Get the list of phases in the chemical system
    auto phases() const -> const PhaseList&;

    /// Calculate the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    auto gibbs_energies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar enthalpies of the species (in units of J/mol).
    auto enthalpies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    auto helmholtz_energies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar entropies of the species (in units of J/K).
    auto entropies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar volumes of the species (in units of m3/mol).
    auto volumes(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar internal energies of the species (in units of J/mol).
    auto internal_energies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto heat_capacities_cp(double T, double P) const -> ThermoVector;

    /// Calculate the concentrations of the species (no uniform units).
    auto concentrations(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activity coefficients of the species.
    auto ln_activity_coefficients(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activities of the species.
    auto ln_activities(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the chemical potentials of the species (in units of J/mol).
    auto chemical_potentials(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the densities of the phases (in units of kg/m3).
    auto densities(double T, double P, const Vector& n) const -> ChemicalVector;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Assemble the formula matrix of a chemical system.
/// The formula matrix of a chemical system is defined as the matrix whose entry
/// `(j, i)` is given by the number of atoms of its `j`-th element in its `i`-th species.
/// @param system The chemical system instance
auto formulaMatrix(const ChemicalSystem& system) -> Matrix;

/// Assemble the balance matrix of a chemical system.
/// The balance matrix of a chemical system is defined as the matrix whose entry
/// `(j, i)` is given by the number of atoms of its `j`-th element in its `i`-th species.
/// The last row of the balance matrix, however, corresponds to the vector of charges of the species.
/// @param system The chemical system instance
auto balanceMatrix(const ChemicalSystem& system) -> Matrix;

} // namespace Reaktor
