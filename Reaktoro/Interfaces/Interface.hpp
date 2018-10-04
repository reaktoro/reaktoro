// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// C++ includes
#include <string>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Thermodynamics/Models/ChemicalModel.hpp>
#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;

/// A class used to interface other codes with Reaktoro.
class Interface
{
public:
    /// Virtual destructor
    virtual ~Interface() = 0;

    /// Return the temperature (in units of K)
    virtual auto temperature() const -> double = 0;

    /// Return the pressure (in units of Pa)
    virtual auto pressure() const -> double = 0;

    /// Return the amounts of the species (in units of mol)
    virtual auto speciesAmounts() const -> Vector = 0;

    /// Return the number of elements
    virtual auto numElements() const -> unsigned = 0;

    /// Return the number of species
    virtual auto numSpecies() const -> unsigned = 0;

    /// Return the number of phases
    virtual auto numPhases() const -> unsigned = 0;

    /// Return the number of species in a phase
    virtual auto numSpeciesInPhase(Index iphase) const -> unsigned = 0;

    /// Return the name of an element
    virtual auto elementName(Index ielement) const -> std::string = 0;

    /// Return the molar mass of an element (in units of kg/mol)
    virtual auto elementMolarMass(Index ielement) const -> double = 0;

    /// Return the stoichiometry of an element in a species
    virtual auto elementStoichiometry(Index ispecies, Index ielement) const -> double = 0;

    /// Return the name of a species
    virtual auto speciesName(Index ispecies) const -> std::string = 0;

    /// Return the name of a phase
    virtual auto phaseName(Index iphase) const -> std::string = 0;

    /// Return the thermodynamic properties of the phases and its species.
    /// @param iphase The index of the phase
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    virtual auto properties(ThermoModelResult& res, double T, double P) -> void = 0;

    /// Return the chemical properties of the phases and its species.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The amounts of the species (in units of mol)
    virtual auto properties(ChemicalModelResult& res, double T, double P, VectorConstRef n) -> void = 0;

    /// Return a clone of this Interface instance.
    virtual auto clone() const -> std::shared_ptr<Interface> = 0;

    /// Return the formula matrix of the species
    auto formulaMatrix() const -> Matrix;

    /// Return the index of an element
    auto indexElement(std::string element) const -> Index;

    /// Return the index of a species
    auto indexSpecies(std::string species) const -> Index;

    /// Return the index of a phase
    auto indexPhase(std::string phase) const -> Index;

    /// Return the index of the phase with a species
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

    /// Return the index of the first species in a phase
    auto indexFirstSpeciesInPhase(Index iphase) const -> Index;

    /// Return a ChemicalSystem instance created from an instance of a class derived from Interface.
    auto system() const -> ChemicalSystem;

    /// Return a ChemicalState instance created from an instance of a class derived from Interface.
    /// @param system The chemical system created using method @ref system.
    auto state(const ChemicalSystem& system) const -> ChemicalState;

    /// Convert the classes derived from Interface into a ChemicalSystem instance
    operator ChemicalSystem() const;
};

} // namespace Reaktoro
