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

    /// Return the residual molar Gibbs energies of the phases
    virtual auto phaseResidualMolarGibbsEnergies() const -> Vector;

    /// Return the residual molar enthalpies of the phases
    virtual auto phaseResidualMolarEnthalpies() const -> Vector;

    /// Return the residual molar isobaric heat capacities of the phases
    virtual auto phaseResidualMolarHeatCapacitiesConstP() const -> Vector;

    /// Return the residual molar isochoric heat capacities of the phases
    virtual auto phaseResidualMolarHeatCapacitiesConstV() const -> Vector;

    /// Return a reference to the low-level Phreeqc instance
    auto phreeqc() -> Phreeqc&;

    /// Return a const reference to the low-level Phreeqc instance
    auto phreeqc() const -> const Phreeqc&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a Phreeqx instance
auto operator<<(std::ostream& out, const Phreeqx& phreeqx) -> std::ostream&;

} // namespace Reaktoro
