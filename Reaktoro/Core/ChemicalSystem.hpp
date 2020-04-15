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

// C++ includes
#include <memory>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ThermoProps;

/// A class to represent a system and its attributes and properties.
/// @see Species, Phase
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a default ChemicalSystem instance
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with given phases.
    explicit ChemicalSystem(const std::vector<Phase>& phases);

    /// Construct a ChemicalSystem instance with given phases and thermodynamic and chemical models.
    // ChemicalSystem(const std::vector<Phase>& phases, const ThermoModel& thermo_model, const ChemicalModel& chemical_model);

    /// Return the number of elements in the system
    auto numElements() const -> unsigned;

    /// Return the number of species in the system
    auto numSpecies() const -> unsigned;

    /// Return the number of species in a phase of the system
    /// @param iphase The index of the phase
    auto numSpeciesInPhase(Index iphase) const -> unsigned;

    /// Return the number of phases in the system
    auto numPhases() const -> unsigned;

    /// Return the list of elements in the system
    auto elements() const -> const std::vector<Element>&;

    /// Return the list of species in the system
    auto species() const -> const std::vector<Species>&;

    /// Return the list of phases in the system
    auto phases() const -> const std::vector<Phase>&;

    // /// Return the thermodynamic model of the system.
    // auto thermoModel() const -> const ThermoModel&;

    // /// Return the chemical model of the system.
    // auto chemicalModel() const -> const ChemicalModel&;

    /// Return the formula matrix of the system
    /// The formula matrix is defined as the matrix whose entry *(j, i)* is
    /// given by the coefficient of the *j*th element in the *i*th species.
    auto formulaMatrix() const -> MatrixXdConstRef;

    /// Return an element of the system
    /// @param index The index of the element
    auto element(Index index) const -> const Element&;

    /// Return an element of the system
    /// @param name The name of the element
    auto element(std::string name) const -> const Element&;

    /// Return a species of the system
    /// @param index The index of the species
    auto species(Index index) const -> const Species&;

    /// Return a species of the system
    /// @param name The name of the species
    auto species(std::string name) const -> const Species&;

    /// Return a phase of the system
    /// @param index The index of the phase
    auto phase(Index index) const -> const Phase&;

    /// Return a phase of the system
    /// @param name The name of the phase
    auto phase(std::string name) const -> const Phase&;

    /// Return the index of an element in the system
    /// @param name The name of the element
    auto indexElement(std::string name) const -> Index;

    /// Return the index of an element in the system. system
    /// It throws an exception if the element does not exist
    /// @param name The name of the element
    auto indexElementWithError(std::string name) const -> Index;

    /// Return the index of a species in the system
    /// @param name The name of the species
    auto indexSpecies(std::string name) const -> Index;

    /// Return the index of a species in the system.
    /// @param name The name of the species
    /// @return The index of the species if found, or a runtime exception otherwise.
    auto indexSpeciesWithError(std::string name) const -> Index;

    /// Return the index of the first species in the system with any of the given names.
    /// @param names The tentative names of the species in the system.
    /// @return The index of the species if found, or the number of species otherwise.
    auto indexSpeciesAny(const std::vector<std::string>& names) const -> Index;

    /// Return the index of the first species in the system with any of the given names.
    /// @param names The tentative names of the species in the system.
    /// @return The index of the species if found, or a runtime exception otherwise.
    auto indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index;

    /// Return the index of a phase in the system
    /// @param name The name of the phase
    auto indexPhase(std::string name) const -> Index;

    /// Return the index of a phase in the system. system
    /// It throws an exception if the phase does not exist
    /// @param name The name of the phase
    auto indexPhaseWithError(std::string name) const -> Index;

    /// Return the index of the phase that contains a given species
    /// @param index The index of the species
    auto indexPhaseWithSpecies(Index index) const -> Index;

    /// Return the index of the first species in a phase
    /// @param iphase The index of the phase
    auto indexFirstSpeciesInPhase(Index iphase) const -> unsigned;

    /// Return the indices of a set of elements in the system
    /// @param names The names of the elements
    auto indicesElements(const std::vector<std::string>& names) const -> Indices;

    /// Return the indices of the elements that compose a species
    /// @param index The index of the species
    auto indicesElementsInSpecies(Index index) const -> Indices;

    /// Return the indices of the elements that compose a set of species
    /// @param indices The indices of the species
    auto indicesElementsInSpecies(const Indices& indices) const -> Indices;

    /// Return the indices of a set of species in the system
    /// @param names The names of the species
    auto indicesSpecies(const std::vector<std::string>& names) const -> Indices;

    /// Return the indices of the species in a given set of phases
    /// @param indices The indices of the phases
    auto indicesSpeciesInPhases(const Indices& indices) const -> Indices;

    /// Return the indices of a set of phases in the system
    /// @param names The names of the phases
    auto indicesPhases(const std::vector<std::string>& names) const -> Indices;

    /// Return the index of the phase that contains a given species
    /// @param indices The indices of the species
    auto indicesPhasesWithSpecies(const Indices& indices) const -> Indices;

    /// Return the indices of the fluid phases.
    auto indicesFluidPhases() const -> Indices;

    /// Return the indices of the species in the fluid phases.
    auto indicesFluidSpecies() const -> Indices;

    /// Return the indices of the solid phases.
    auto indicesSolidPhases() const -> Indices;

    /// Return the indices of the species in the solid phases.
    auto indicesSolidSpecies() const -> Indices;

    /// Calculate the standard thermodynamic properties of the species.
    /// @param T The temperature of the system (in K)
    /// @param P The pressure of the system (in Pa)
    auto props(double T, double P) const -> ThermoProps;

    /// Calculate the thermodynamic and chemical properties of the chemical system.
    /// @param T The temperature of the system (in K)
    /// @param P The pressure of the system (in Pa)
    /// @param n The amounts of the species in the system (in mol)
    auto props(double T, double P, VectorXrConstRef n) const -> ChemicalProps;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a ChemicalSystem instance
auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&;

} // namespace Reaktoro
