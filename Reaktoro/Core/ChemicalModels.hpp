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
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

class ChemicalModels
{
public:
    /// Construct a default ChemicalModels instance
    ChemicalModels();

    /// Construct a copy of a ChemicalModels instance
    ChemicalModels(const ChemicalModels& other);

    /// Destroy a ChemicalModels instance
    virtual ~ChemicalModels();

    /// Assign a ChemicalModels instance to this
    auto operator=(ChemicalModels other) -> ChemicalModels&;

    /// Set the function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    auto setStandardGibbsEnergyFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the apparent standard molar enthalpies of the species (in units of J/mol).
    auto setStandardEnthalpyFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    auto setStandardHelmholtzEnergyFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the apparent standard molar internal energies of the species (in units of J/mol).
    auto setStandardInternalEnergyFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the standard molar entropies of the species (in units of J/K).
    auto setStandardEntropyFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the standard molar volumes of the species (in units of m3/mol).
    auto setStandardVolumeFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto setStandardHeatCapacityFunction(const ThermoVectorFunction& function) -> void;

    /// Set the function for the concentrations of the species (no uniform units).
    auto setConcentrationFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the natural log of the activity coefficients of the species.
    auto setActivityCoefficientFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the natural log of the activities of the species.
    auto setActivityFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the chemical potentials of the species (in units of J/mol).
    auto setChemicalPotentialFunction(const ChemicalVectorFunction& function) -> void;

    /// Set the function for the molar volumes of the phases (in units of m3/mol).
    auto setPhaseMolarVolumeFunction(const ChemicalVectorFunction& function) -> void;

    /// Return the function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    auto standardGibbsEnergyFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the apparent standard molar enthalpies of the species (in units of J/mol).
    auto standardEnthalpyFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    auto standardHelmholtzEnergyFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the apparent standard molar internal energies of the species (in units of J/mol).
    auto standardInternalEnergyFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the standard molar entropies of the species (in units of J/K).
    auto standardEntropyFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the standard molar volumes of the species (in units of m3/mol).
    auto standardVolumeFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    auto standardHeatCapacityFunction() const -> const ThermoVectorFunction&;

    /// Return the function for the concentrations of the species (no uniform units).
    auto concentrationFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the natural log of the activity coefficients of the species.
    auto activityCoefficientFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the natural log of the activities of the species.
    auto activityFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the chemical potentials of the species (in units of J/mol).
    auto chemicalPotentialFunction() const -> const ChemicalVectorFunction&;

    /// Return the function for the molar volumes of the phases (in units of m3/mol).
    auto phaseMolarVolumeFunction() const -> const ChemicalVectorFunction&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
