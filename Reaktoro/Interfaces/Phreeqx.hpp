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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

// Forward declarations
class Phreeqc;

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class ChemicalState;

class Phreeqx
{
public:
    /// Construct a default Phreeqx instance
    Phreeqx();

    /// Construct a Phreeqx instance from a specification file
    /// @param database The name of the database file
    /// @param script The name of the script file
    Phreeqx(std::string database, std::string script);

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

    /// Return the standard molar Gibbs free energy of a species (in units of J/mol)
    virtual auto standardMolarGibbsEnergy(Index ispecies) const -> double;

    /// Return the standard molar enthalpy of a species (in units of J/mol)
    virtual auto standardMolarEnthalpy(Index ispecies) const -> double;

    /// Return the standard molar volume of a species (in units of m3/mol)
    virtual auto standardMolarVolume(Index ispecies) const -> double;

    /// Return the standard molar isobaric heat capacity of a species (in units of J/(mol*K))
    virtual auto standardMolarHeatCapacityConstP(Index ispecies) const -> double;

    /// Return the standard molar isochoric heat capacity of a species (in units of J/(mol*K))
    virtual auto standardMolarHeatCapacityConstV(Index ispecies) const -> double;

    /// Return the ln activity coefficients of the species in a phase.
    virtual auto lnActivityCoefficients(Index iphase) const -> Vector;

    /// Return the ln activities of the species in a phase
    virtual auto lnActivities(Index iphase) const -> Vector;

    /// Return the molar volume of a phase
    virtual auto phaseMolarVolume(Index iphase) const -> double;

    /// Return the residual molar Gibbs energy of a phase
    virtual auto phaseResidualMolarGibbsEnergy(Index iphase) const -> double;

    /// Return the residual molar enthalpy of a phase
    virtual auto phaseResidualMolarEnthalpy(Index iphase) const -> double;

    /// Return the residual molar isobaric heat capacity of a phase
    virtual auto phaseResidualMolarHeatCapacityConstP(Index iphase) const -> double;

    /// Return the residual molar isochoric heat capacity of a phase
    virtual auto phaseResidualMolarHeatCapacityConstV(Index iphase) const -> double;


















    /// Set the temperature of the Phreeqx instance (in units of K)
    auto setTemperature(double val) -> void;

    /// Set the pressure of the Phreeqx instance (in units of Pa)
    auto setPressure(double val) -> void;

    /// Set the amounts of the species of the Phreeqx instance (in units of mol)
    auto setSpeciesAmounts(const Vector& n) -> void;

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
    auto indexPhaseWithSpecies(unsigned ispecies) const -> unsigned;

    /// Return the number of atoms of an element in a species
    /// @param ielement The index of the element
    /// @param ispecies The index of the species
    auto elementCoefficientInSpecies(unsigned ielement, unsigned ispecies) const -> double;

    /// Return the indices and number of atoms of the elements that compose a species
    /// @param index The index of the species
    auto elementsInSpecies(unsigned index) const -> std::map<unsigned, double>;

    /// Return the molar mass of an element (in units of kg/mol)
    /// @param index The index of the element
    auto elementMolarMass(unsigned index) const -> double;

    /// Return the temperature of the Phreeqx instance (in units of K)
    auto temperature() const -> double;

    /// Return the pressure of the Phreeqx instance (in units of Pa)
    auto pressure() const -> double;

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

    /// Return the standard molar Gibbs free energies of the species (in units of J/mol)
    auto standardGibbsEnergies() -> Vector;

    /// Return the standard molar volumes of the species (in units of m3/mol)
    auto standardVolumes() -> Vector;

    /// Return the natural logarithm of the activities of the species
    auto activities() -> Vector;

    /// Return the chemical potentials of the species (in units of J/mol)
    auto chemicalPotentials() -> Vector;

    /// Return the molar volumes of the phases (in units of m3/mol)
    auto phaseMolarVolumes() -> Vector;

    /// Return a reference to the low-level Phreeqc instance
    auto phreeqc() -> Phreeqc&;

    /// Return a const reference to the low-level Phreeqc instance
    auto phreeqc() const -> const Phreeqc&;

    /// Convert this Phreeqx instance into a ChemicalSystem instance
    operator ChemicalSystem() const;

    /// Convert this Phreeqx instance into a ChemicalState instance
    operator ChemicalState() const;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a Phreeqx instance
auto operator<<(std::ostream& out, const Phreeqx& phreeqx) -> std::ostream&;

} // namespace Reaktoro
