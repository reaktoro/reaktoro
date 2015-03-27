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
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Phase.hpp>

namespace Reaktor {

/// A struct to represent a more detailed model configuration of a Multiphase object.
/// @see Multiphase, Phase
/// @ingroup Core
struct MultiphaseModel
{
    /// The function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_gibbs_energy_fn;

    /// The function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_helmholtz_energy_fn;

    /// The function for the apparent standard molar internal energies of the species (in units of J/mol).
    ThermoVectorFunction standard_internal_energy_fn;

    /// The function for the apparent standard molar enthalpies of the species (in units of J/mol).
    ThermoVectorFunction standard_enthalpy_fn;

    /// The function for the standard molar entropies of the species (in units of J/K).
    ThermoVectorFunction standard_entropy_fn;

    /// The function for the standard molar volumes of the species (in units of m3/mol).
    ThermoVectorFunction standard_volume_fn;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoVectorFunction standard_heat_capacity_fn;

    /// The function for the concentrations of the species (no uniform units).
    ChemicalVectorFunction concentration_fn;

    /// The function for the natural log of the activity coefficients of the species.
    ChemicalVectorFunction activity_coefficient_fn;

    /// The function for the natural log of the activities of the species.
    ChemicalVectorFunction activity_fn;

    /// The function for the chemical potentials of the species (in units of J/mol).
    ChemicalVectorFunction chemical_potential_fn;

    /// The function for the molar volumes of the phases (in units of m3/mol).
    ChemicalVectorFunction phase_molar_volume_fn;
};

/// A class to represent a multiphase system and its attributes and properties.
/// @see Species, Phase
/// @ingroup Core
class Multiphase
{
public:
    /// Construct a default Multiphase instance
    Multiphase();

    /// Construct a Multiphase instance with given phases.
    explicit Multiphase(const std::vector<Phase>& phases);

    /// Construct a Multiphase instance with given phases and model configuration.
    explicit Multiphase(const std::vector<Phase>& phases, const MultiphaseModel& model);

    /// Destroy this Multiphase instance
    virtual ~Multiphase();

    /// Return the number of elements in the multiphase system
    auto numElements() const -> unsigned;

    /// Return the number of species in the multiphase system
    auto numSpecies() const -> unsigned;

    /// Return the number of species in a phase of the multiphase system
    /// @param iphase The index of the phase
    auto numSpeciesInPhase(Index iphase) const -> unsigned;

    /// Return the number of phases in the multiphase system
    auto numPhases() const -> unsigned;

    /// Return the list of elements in the multiphase system
    auto elements() const -> const std::vector<Element>&;

    /// Return the list of species in the multiphase system
    auto species() const -> const std::vector<Species>&;

    /// Return the list of phases in the multiphase system
    auto phases() const -> const std::vector<Phase>&;

    /// Return the formula matrix of the multiphase system
    /// The formula matrix is defined as the matrix whose entry `(j, i)`
    /// is given by the number of atoms of its `j`-th element in its `i`-th species.
    auto formulaMatrix() const -> const Matrix&;

    /// Return an element of the multiphase system
    /// @param index The index of the element
    auto element(Index index) const -> const Element&;

    /// Return an element of the multiphase system
    /// @param name The name of the element
    auto element(std::string name) const -> const Element&;

    /// Return a species of the multiphase system
    /// @param index The index of the species
    auto species(Index index) const -> const Species&;

    /// Return a species of the multiphase system
    /// @param name The name of the species
    auto species(std::string name) const -> const Species&;

    /// Return a phase of the multiphase system
    /// @param index The index of the phase
    auto phase(Index index) const -> const Phase&;

    /// Return a phase of the multiphase system
    /// @param name The name of the phase
    auto phase(std::string name) const -> const Phase&;

    /// Return the index of an element in the multiphase system
    /// @param name The name of the element
    auto indexElement(std::string name) const -> Index;

    /// Return the index of an element in the multiphase. system
    /// It throws an exception if the element does not exist
    /// @param name The name of the element
    auto indexElementWithError(std::string name) const -> Index;

    /// Return the index of a species in the multiphase system
    /// @param name The name of the species
    auto indexSpecies(std::string name) const -> Index;

    /// Return the index of a species in the multiphase. system
    /// It throws an exception if the species does not exist
    /// @param name The name of the species
    auto indexSpeciesWithError(std::string name) const -> Index;

    /// Return the index of a phase in the multiphase system
    /// @param name The name of the phase
    auto indexPhase(std::string name) const -> Index;

    /// Return the index of a phase in the multiphase. system
    /// It throws an exception if the phase does not exist
    /// @param name The name of the phase
    auto indexPhaseWithError(std::string name) const -> Index;

    /// Return the index of the phase that contains a given species
    /// @param index The index of the species
    auto indexPhaseWithSpecies(Index index) const -> Index;

    /// Return the indices of a set of elements in the multiphase system
    /// @param name The names of the elements
    auto indicesElements(const std::vector<std::string>& names) const -> Indices;

    /// Return the indices of a set of species in the multiphase system
    /// @param names The names of the species
    auto indicesSpecies(const std::vector<std::string>& names) const -> Indices;

    /// Return the indices of a set of phases in the multiphase system
    /// @param names The names of the phases
    auto indicesPhases(const std::vector<std::string>& names) const -> Indices;

    /// Return the index of the phase that contains a given species
    /// @param indices The indices of the species
    auto indicesPhasesWithSpecies(const Indices& indices) const -> Indices;

    /// Return the indices of the elements that compose a species
    /// @param index The index of the species
    auto indicesElementsInSpecies(Index index) const -> Indices;

    /// Return the indices of the elements that compose a set of species
    /// @param indices The indices of the species
    auto indicesElementsInSpecies(const Indices& indices) const -> Indices;

    /// Return the index of the first species in a phase
    /// @param The index of the phase
    auto indexFirstSpeciesInPhase(Index iphase) const -> unsigned;

    /// Calculate the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardGibbsEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar enthalpies of the species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardEnthalpies(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardHelmholtzEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar entropies of the species (in units of J/K).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardEntropies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar volumes of the species (in units of m3/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardVolumes(double T, double P) const -> ThermoVector;

    /// Calculate the apparent standard molar internal energies of the species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardInternalEnergies(double T, double P) const -> ThermoVector;

    /// Calculate the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto standardHeatCapacities(double T, double P) const -> ThermoVector;

    /// Calculate the concentrations of the species (no uniform units).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto concentrations(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activity coefficients of the species.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param n The molar amounts of the species (in units of mol)
    auto activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the natural log of the activities of the species.
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param n The molar amounts of the species (in units of mol)
    auto activities(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the chemical potentials of the species (in units of J/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param n The molar amounts of the species (in units of mol)
    auto chemicalPotentials(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the molar volumes of the phases (in units of m3/mol).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    /// @param n The molar amounts of the species (in units of mol)
    auto phaseMolarVolumes(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the densities of the phases (in units of kg/m3).
    /// @param T The temperature value (in units of K)
    /// @param P The pressure value (in units of Pa)
    auto phaseDensities(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Return the total molar amounts in each phase (in units of mol)
    /// @param n The molar amounts of the species (in units of mol)
    auto phaseMolarAmounts(const Vector& n) const -> ChemicalVector;

    /// Return the total mass amounts in each phase (in units of kg)
    /// @param n The molar amounts of the species (in units of mol)
    auto phaseMassAmounts(const Vector& n) const -> Vector;

    /// Calculate the volumes of the phases (in units of m3).
    /// @param n The molar amounts of the species (in units of mol)
    auto phaseVolumes(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the molar amounts of the elements (in units of mol)
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmounts(const Vector& n) const -> Vector;

    /// Calculate the molar amounts of the elements in a given phase (in units of mol)
    /// @param iphase The index of the phase
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmountsInPhase(Index iphase, const Vector& n) const -> Vector;

    /// Calculate the molar amounts of the elements in a given set of species (in units of mol)
    /// @param ispecies The indices of the species
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmountsInSpecies(const Indices& ispecies, const Vector& n) const -> Vector;

    /// Calculate the molar amount of an elements (in units of mol)
    /// @param ielement The index of the element
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmount(Index ielement, const Vector& n) const -> double;

    /// Calculate the molar amounts of the elements in a given phase (in units of mol)
    /// @param ielement The index of the element
    /// @param iphase The index of the phase
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmountInPhase(Index ielement, Index iphase, const Vector& n) const -> double;

    /// Calculate the molar amounts of the elements in a given set of species (in units of mol)
    /// @param ielement The index of the element
    /// @param ispecies The indices of the species in the set
    /// @param n The molar amounts of the species (in units of mol)
    auto elementAmountInSpecies(Index ielement, const Indices& ispecies, const Vector& n) const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a Multiphase instance
auto operator<<(std::ostream& out, const Multiphase& multiphase) -> std::ostream&;

} // namespace Reaktor
