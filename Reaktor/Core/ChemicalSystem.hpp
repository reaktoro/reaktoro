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
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Component.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Phase.hpp>

namespace Reaktor {

/// The type used to define the model functions of a ChemicalSystem instance
/// @see ChemicalSystem
/// @ingroup Core
struct ChemicalSystemModels
{
    /// The thermodynamic model function for the calculation of the
    /// apparent standard molar Gibbs free energies of the species (in units of J/mol).
    ThermoVectorFunction g0;

    /// The thermodynamic model function for the calculation of the
    /// apparent standard molar enthalpies of the species (in units of J/mol).
    ThermoVectorFunction h0;

    /// The thermodynamic model function for the calculation of the
    /// apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    ThermoVectorFunction a0;

    /// The thermodynamic model function for the calculation of the
    /// standard molar entropies of the species (in units of J/K).
    ThermoVectorFunction s0;

    /// The thermodynamic model function for the calculation of the
    /// standard molar volumes of the species (in units of m3/mol).
    ThermoVectorFunction v0;

    /// The thermodynamic model function for the calculation of the
    /// apparent standard molar internal energies of the species (in units of J/mol).
    ThermoVectorFunction u0;

    /// The thermodynamic model function for the calculation of the
    /// the standard molar isobaric heat capacity of the species (in units of J/(mol*K))
    ThermoVectorFunction cp0;

    /// The chemical model function for the calculation of the
    /// concentrations of the species (no uniform units)
    ChemicalVectorFunction c;

    /// The chemical model function for the calculation of the
    /// natural log of the activity coefficients of the species (in units of J/mol)
    ChemicalVectorFunction ln_gamma;

    /// The chemical model function for the calculation of the
    /// natural log of the activities of the species (in units of J/mol)
    ChemicalVectorFunction ln_a;

    /// The chemical model function for the calculation of the
    /// molar Gibbs energies of the species (in units of J/mol)
    ChemicalVectorFunction g;

    /// The chemical model function for the calculation of the
    /// molar densities of the phases (in units of mol/kg3)
    ChemicalVectorFunction xi;

    /// The chemical model function for the calculation of the
    /// densities of the phases (in units of kg/m3)
    ChemicalVectorFunction rho;
};

/// The type used to define the attributes and model functions of a ChemicalSystem instance
/// @see ChemicalSystem
/// @ingroup Core
struct ChemicalSystemData
{
    /// The list of components in the chemical system
    ComponentList components;

    /// The list of species in the chemical system
    SpeciesList species;

    /// The list of phases in the chemical system
    PhaseList phases;

    /// The thermodynamic and chemical models of the species and phases in the chemical system
    ChemicalSystemModels models;
};

/// The type used to define a chemical system and its attributes
/// @see Species, Phase
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a default ChemicalSystem instance
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with all its attributes
    ChemicalSystem(const ChemicalSystemData& data);

    /// Get the list of components in the chemical system
    auto components() const -> const ComponentList&;

    /// Get the list of species in the chemical system
    auto species() const -> const SpeciesList&;

    /// Get the list of phases in the chemical system
    auto phases() const -> const PhaseList&;

    /// Get the thermodynamic and chemical models of the chemical system
    auto models() const -> const ChemicalSystemModels&;

private:
    /// The immutable shared data of the ChemicalSystem class
    std::shared_ptr<ChemicalSystemData> data;
};

} // namespace Reaktor
