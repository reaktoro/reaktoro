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

namespace Reaktoro {

class Interface
{
public:
    Interface();

    virtual ~Interface() = 0;

    /// Return the number of elements
    auto numElements() const -> unsigned = 0;

    /// Return the number of species
    auto numSpecies() const -> unsigned = 0;

    /// Return the number of phases
    auto numPhases() const -> unsigned = 0;

    /// Return the number of species in a phase
    auto numSpeciesInPhase(Index iphase) const -> unsigned = 0;

    /// Return the name of an element
    auto elementName(Index ielement) const -> std::string = 0;

    /// Return the molar mass of an element (in units of kg/mol)
    auto elementMolarMass(Index ielement) const -> double = 0;

    /// Return the name of a species
    auto speciesName(Index ispecies) const -> std::string = 0;

    /// Return the molar mass of a species (in units of kg/mol)
    auto speciesMolarMass(Index ispecies) const -> double = 0;

    /// Return the electrical charge of a species
    auto speciesCharge(Index ispecies) const -> double = 0;

    /// Return the indices and stoichiometries of the elements in a species
    auto speciesElements(Index ispecies) const -> std::map<Index, double> = 0;

    /// Return the name of a phase
    auto phaseName(Index iphase) const -> std::string = 0;




    /// Return the thermodynamic state of a phase.
    auto phaseThermoState(Index iphase, double T, double P) -> PhaseThermoModelResult = 0;

    /// Return the chemical state of a phase.
    auto phaseChemicalState(Index iphase, double T, double P, const Vector& n) -> PhaseChemicalModelResult = 0;
};

} // namespace Reaktoro
