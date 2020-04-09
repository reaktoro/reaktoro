// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Interfaces/Interface.hpp>

// Forward declarations
class TNode;

namespace Reaktoro {

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

    /// Return the temperature (in units of K)
    virtual auto temperature() const -> double;

    /// Return the pressure (in units of Pa)
    virtual auto pressure() const -> double;

    /// Return the amounts of the species (in units of mol)
    virtual auto speciesAmounts() const -> VectorXr;

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
    virtual auto properties(ChemicalModelResult& res, double T, double P, VectorXrConstRef n) -> void;

    /// Return a clone of this Gems instance
    virtual auto clone() const -> std::shared_ptr<Interface>;

    /// Set the temperature and pressure of the Gems instance.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    auto set(double T, double P) -> void;

    /// Set the temperature, pressure and molar composition of the Gems instance.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The composition of the species (in units of mol)
    auto set(double T, double P, VectorXrConstRef n) -> void;

    /// Set the options of the Gems instance
    auto setOptions(const GemsOptions& options) -> void;

    /// Calculate the equilibrium state of the system
    /// @param T The temperature for the equilibrium calculation (in units of K)
    /// @param P The pressure for the equilibrium calculation (in units of Pa)
    /// @param n The amounts of the elements (in units of mol)
    auto equilibrate(double T, double P, VectorXrConstRef b) -> void;

    /// Return the convergence result of the equilibrium calculation
    auto converged() const -> bool;

    /// Return the number of iterations of the equilibrium calculation
    auto numIterations() const -> unsigned;

    /// Return the wall time of the equilibrium calculation (in units of s)
    auto elapsedTime() const -> double;

    /// Return a shared pointer to the TNode instance of Gems
    auto node() const -> std::shared_ptr<TNode>;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
