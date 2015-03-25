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

// C++ includes
#include <map>
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>

// Forward declarations
class TNode;

namespace Reaktor {

// Forward declarations
class ChemicalSystem;
class ChemicalState;

/// A type that describes the options for Gems
struct GemsOptions
{
    /// The flag that indicates if smart start initial approximation is used
    bool warmstart = true;
};

/// A wrapper class for Gems code
class Gems
{
public:
    /// Construct a default Gems instance
    Gems();

    /// Construct a Gems instance from a specification file
    /// @param filename The name of the file containing the definition of the chemical system
    Gems(std::string filename);

    /// Set the temperature of the Gems instance (in units of K)
    auto setTemperature(double val) -> void;

    /// Set the pressure of the Gems instance (in units of Pa)
    auto setPressure(double val) -> void;

    /// Set the amounts of the species of the Gems instance (in units of mol)
    auto setSpeciesAmounts(const Vector& n) -> void;

    /// Set the amounts of the elements of the Gems instance (in units of mol)
    auto setElementAmounts(const Vector& b) -> void;

    /// Set the options of the Gems instance
    auto setOptions(const GemsOptions& options) -> void;

    /// Return the number of elements
    auto numElements() const -> unsigned;

    /// Return the number of species
    auto numSpecies() const -> unsigned;

    /// Return the number of phases
    auto numPhases() const -> unsigned;

    /// Return the number of species in a phase
    /// @param index The index of the phase
    auto numSpeciesInPhase(unsigned index) const -> unsigned;

    /// Return the name of an element
    /// @param index The index of the element
    auto elementName(unsigned index) const -> std::string;

    /// Return the name of a species
    /// @param index The index of the species
    auto speciesName(unsigned index) const -> std::string;

    /// Return the name of a phase
    /// @param index The index of the phase
    auto phaseName(unsigned index) const -> std::string;

    /// Return the index of an element
    /// @param name The name of the element
    auto indexElement(std::string name) const -> unsigned;

    /// Return the index of a species
    /// @param name The name of the species
    auto indexSpecies(std::string name) const -> unsigned;

    /// Return the index of a phase
    /// @param name The name of the phase
    auto indexPhase(std::string name) const -> unsigned;

    /// Return the index of the phase with a species
    /// @param ispecies The index of the species
    auto indexPhaseWithSpecies(unsigned ispecies) const -> Index;

    /// Return the number of atoms of an element in a species
    /// @param ielement The index of the element
    /// @param ispecies The index of the species
    auto elementAtomsInSpecies(unsigned ielement, unsigned ispecies) const -> double;

    /// Return the electrical charge of a species
    /// @param index The index of the species
    auto speciesCharge(unsigned index) const -> double;

    /// Return the indices and number of atoms of the elements that compose a species
    /// @param index The index of the species
    auto elementsInSpecies(unsigned index) const -> std::map<unsigned, double>;

    /// Return the molar mass of an element (in units of kg/mol)
    /// @param index The index of the element
    auto elementMolarMass(unsigned index) const -> double;

    /// Return the molar mass of a species (in units of kg/mol)
    /// @param index The index of the species
    auto speciesMolarMass(unsigned index) const -> double;

    /// Return the temperature of the Gems instance (in units of K)
    auto temperature() const -> double;

    /// Return the pressure of the Gems instance (in units of Pa)
    auto pressure() const -> double;

    /// Return the amounts of the elements (in units of mol)
    auto elementAmounts() const -> Vector;

    /// Return the amounts of the species (in units of mol)
    auto speciesAmounts() const -> Vector;

    /// Return the amounts of a species (in units of mol)
    /// @param index The index of the species
    auto speciesAmount(unsigned index) const -> double;

    /// Return the amounts of the species in a given phase (in units of mol)
    /// @param index The index of the phase
    auto speciesAmountsInPhase(unsigned index) const -> Vector;

    /// Return the formula matrix of the species
    auto formulaMatrix() const -> Matrix;

    /// Return the molar standard Gibbs free energies of the species
    auto standardGibbsEnergies() -> Vector;

    /// Return the standard molar volumes of the species (in units of m3/mol)
    auto standardVolumes() -> Vector;

    /// Return the chemical potentials of the species
    auto chemicalPotentials() -> Vector;

    /// Return the molar volumes of the phases (in units of m3/mol)
    auto phaseMolarVolumes() -> Vector;

    /// Calculate the equilibrium state of the system
    auto equilibrate() -> void;

    /// Return the convergence result of the equilibrium calculation
    auto converged() const -> bool;

    /// Return the number of iterations of the equilibrium calculation
    auto numIterations() const -> unsigned;

    /// Return the wall time of the equilibrium calculation (in units of s)
    auto elapsedTime() const -> double;

    /// Return a reference to the TNode instance of Gems
    auto node() -> TNode&;

    /// Return a const reference to the TNode instance of Gems
    auto node() const -> const TNode&;

    /// Convert this Gems instance into a ChemicalSystem instance
    operator ChemicalSystem() const;

    /// Convert this Gems instance into a ChemicalState instance
    operator ChemicalState() const;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktor
