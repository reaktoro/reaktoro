// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Interfaces/Interface.hpp>

// Forward declarations
class PHREEQC;

namespace Reaktoro {

// Forward declarations
class ReactionEquation;

class Phreeqc : public Interface
{
public:
    /// Construct a default Phreeqc instance
    Phreeqc();

    /// Construct a Phreeqc instance with given PHREEQC database.
    /// @param database The name of the database file
    Phreeqc(std::string database);

    /// Destroy this Phreeqc instance
    virtual ~Phreeqc();

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

    /// Return the thermodynamic properties of the phases and its species.
    /// @param iphase The index of the phase
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    virtual auto properties(ThermoModelResult& res, double T, double P) -> void;

    /// Return the chemical properties of the phases and its species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The amounts of the species (in units of mol)
    virtual auto properties(ChemicalModelResult& res, double T, double P, VectorConstRef n) -> void;

    /// Return a clone of this Phreeqc instance
    virtual auto clone() const -> std::shared_ptr<Interface>;

    /// Set the temperature and pressure of the interfaced code.
    /// This method should be used to update all thermodynamic properties
    /// that depend only on temperature and pressure, such as standard thermodynamic
    /// properties of the species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    auto set(double T, double P) -> void;

    /// Set the temperature, pressure and species composition of the interfaced code.
    /// This method should be used to update all thermodynamic properties
    /// that depend only on temperature and pressure, such as standard thermodynamic
    /// properties of the species, as well as chemical properties that depend on the
    /// composition of the species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The composition of the species (in units of mol)
    auto set(double T, double P, VectorConstRef n) -> void;

    /// Load a PHREEQC database.
    /// This method will initialize the Phreeqc instance with all species and reactions
    /// found in the given database.
    /// To execute a Phreeqc script file, check the `execute` method.
    /// @see execute
    /// @param database The path to the database file.
    auto load(std::string database) -> void;

    /// Execute a PHREEQC input script either provided as a file or input string.
    /// This method will execute the given PHREEQC input script and put this Phreeqc
    /// instance in a state with active species and phases from the last PHREEQC
    /// calculation specified in the script file. This can then be used to
    /// initialize a ChemicalSystem instance with such configuration.
    /// @param input The input either as a filename or as an input script coded in a string.
    /// @param output The name of the file where the result should be output
    auto execute(std::string input, std::string output) -> void;

    /// Execute a PHREEQC input script either provided as a file or input string.
    /// @param input The input either as a filename or as an input script coded in a string.
    auto execute(std::string input) -> void;

    /// Reset this Phreeqc instance to a clean state.
    /// After this Phreeqc instance is reset, one must again
    /// load a new database and execute a new input script file.
    /// @see load, execute
    auto reset() -> void;

    /// Return the system of reactions.
    auto reactions() const -> std::vector<ReactionEquation>;

    /// Return the stoichiometric matrix of the system of reactions.
    auto stoichiometricMatrix() const -> Matrix;

    /// Return the standard molar Gibbs free energies of the species (in units of J/mol)
    auto standardMolarGibbsEnergies() const -> Vector;

    /// Return the standard molar enthalpies of the species (in units of J/mol)
    auto standardMolarEnthalpies() const -> Vector;

    /// Return the standard molar volumes of the species (in units of m3/mol)
    auto standardMolarVolumes() const -> Vector;

    /// Return the standard molar isobaric heat capacities of the species (in units of J/(mol*K))
    auto standardMolarHeatCapacitiesConstP() const -> Vector;

    /// Return the standard molar isochoric heat capacities of the species (in units of J/(mol*K))
    auto standardMolarHeatCapacitiesConstV() const -> Vector;

    /// Return the ln activity coefficients of the species
    auto lnActivityCoefficients() const -> Vector;

    /// Return the ln activity contants of the species
    auto lnActivityConstants() const -> Vector;

    /// Return the ln activities of the species
    auto lnActivities() const -> Vector;

    /// Return the ln equilibrium constants of the reactions
    auto lnEquilibriumConstants() const -> Vector;

    /// Return the molar volumes of the phases
    auto phaseMolarVolumes() const -> Vector;

    /// Return a reference to the low-level Phreeqc instance
    auto phreeqc() -> PHREEQC&;

    /// Return a const reference to the low-level Phreeqc instance
    auto phreeqc() const -> const PHREEQC&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
