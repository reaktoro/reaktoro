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
#include <map>
#include <memory>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Interfaces/Interface.hpp>

// Forward declarations
class TNode;

namespace Reaktoro {

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
class Gems : public Interface
{
public:
    /// Construct a default Gems instance
    Gems();

    /// Construct a Gems instance from a specification file
    /// @param filename The name of the file containing the definition of the chemical system
    Gems(std::string filename);

    /// Destroy this Gems instance
    virtual ~Gems();

    /// Set the temperature and pressure of the interfaced code.
    /// This method should be used to update all thermodynamic properties
    /// that depend only on temperature and pressure, such as standard thermodynamic
    /// properties of the species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    virtual auto set(double T, double P) -> void;

    /// Set the temperature, pressure and species composition of the interfaced code.
    /// This method should be used to update all thermodynamic properties
    /// that depend only on temperature and pressure, such as standard thermodynamic
    /// properties of the species, as well as chemical properties that depend on the
    /// composition of the species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The composition of the species (in units of mol)
    virtual auto set(double T, double P, const Vector& n) -> void;

    /// Return the temperature (in units of K)
    virtual auto temperature() const -> double;

    /// Return the pressure (in units of Pa)
    virtual auto pressure() const -> double;

    /// Return the amounts of the species (in units of mol)
    virtual auto speciesAmounts() const -> Vector;

    /// Return the number of elements
    virtual auto numElements() const -> unsigned;

    /// Return the number of species
    virtual auto numSpecies() const -> unsigned;

    /// Return the number of phases
    virtual auto numPhases() const -> unsigned;

    /// Return the number of species in a phase
    virtual auto numSpeciesInPhase(Index iphase) const -> unsigned;

    /// Return the name of an element
    virtual auto elementName(Index ielement) const -> std::string;

    /// Return the molar mass of an element (in units of kg/mol)
    virtual auto elementMolarMass(Index ielement) const -> double;

    /// Return the stoichiometry of an element in a species
    virtual auto elementStoichiometry(Index ispecies, Index ielement) const -> double;

    /// Return the name of a species
    virtual auto speciesName(Index ispecies) const -> std::string;

    /// Return the name of a phase
    virtual auto phaseName(Index iphase) const -> std::string;

    /// Return the standard reference state of a phase (`"IdealGas"` or `"IdealSolution`)
    virtual auto phaseReferenceState(Index iphase) const -> std::string;

    /// Return the standard molar Gibbs free energies of the species (in units of J/mol)
    virtual auto standardMolarGibbsEnergies() const -> Vector;

    /// Return the standard molar enthalpies of the species (in units of J/mol)
    virtual auto standardMolarEnthalpies() const -> Vector;

    /// Return the standard molar volumes of the species (in units of m3/mol)
    virtual auto standardMolarVolumes() const -> Vector;

    /// Return the standard molar isobaric heat capacities of the species (in units of J/(mol*K))
    virtual auto standardMolarHeatCapacitiesConstP() const -> Vector;

    /// Return the standard molar isochoric heat capacities of the species (in units of J/(mol*K))
    virtual auto standardMolarHeatCapacitiesConstV() const -> Vector;

    /// Return the ln activity coefficients of the species
    virtual auto lnActivityCoefficients() const -> Vector;

    /// Return the ln activities of the species
    virtual auto lnActivities() const -> Vector;

    /// Return the molar volumes of the phases
    virtual auto phaseMolarVolumes() const -> Vector;

    /// Set the options of the Gems instance
    auto setOptions(const GemsOptions& options) -> void;

    /// Calculate the equilibrium state of the system
    /// @param T The temperature for the equilibrium calculation (in units of K)
    /// @param P The pressure for the equilibrium calculation (in units of Pa)
    /// @param n The amounts of the elements (in units of mol)
    auto equilibrate(double T, double P, const Vector& b) -> void;

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

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
