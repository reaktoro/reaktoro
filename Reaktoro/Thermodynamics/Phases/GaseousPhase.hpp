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
#include <string>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivity.hpp>

namespace Reaktoro {

// Forward declarations
class GaseousMixture;

/// Class that defines a gaseous phase
class GaseousPhase : public Phase
{
public:
    /// Construct a default GaseousPhase instance.
    GaseousPhase();

    /// Construct a copy of a GaseousPhase instance
    GaseousPhase(const GaseousPhase& other);

    /// Construct an GaseousPhase instance with given gaseous mixture.
    explicit GaseousPhase(const GaseousMixture& mixture);

    /// Destroy the GaseousPhase instance.
    virtual ~GaseousPhase();

    /// Assign a GaseousPhase instance to this
    auto operator=(GaseousPhase other) -> GaseousPhase&;

    /// Set the chemical model of the phase with the ideal gas equation of state.
    auto setChemicalModelIdeal() -> void;

    /// Set the chemical model of the phase with the Spyecher et al. (2003) equation of state.
    /// This model only supports the gaseous species `H2O(g)` and `CO2(g)`. Any other species will result in
    /// a runtime error.
    /// Reference: *Spycher, N., Pruess, K., Ennis-King, J. (2003). CO2-H2O mixtures in the
    /// geological sequestration of CO2. I. Assessment and calculation of mutual solubilities from 12 to 100°C
    /// and up to 600 bar. Geochimica et Cosmochimica Acta, 67(16), 3015–3031*.
    auto setChemicalModelSpycherEtAl2003() -> void;

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see GaseousActivityFunction
    auto setActivityModel(std::string species, const GaseousActivityFunction& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(std::string species) -> void;

    /// Set the activity model of CO2(g) to be the one of Duan et al. (2006).
    auto setActivityModelDuanSunCO2() -> void;

    /// Set the activity model of H2O(g) and CO2(g) to be the one of Spycher et al. (2003).
    auto setActivityModelSpycherPruessH2OCO2() -> void;

    /// Set the activity model of H2O(g), CO2(g) and CH4(g) to be the one of Spycher and Reed (1988).
    auto setActivityModelSpycherReedH2OCO2CH4() -> void;

    /// Set the activity model of a gaseous species to be the one of Peng and Robinson (1978).
    /// @param species The name of the gaseous species
    auto setActivityModelPengRobinson(std::string species) -> void;

    /// Return the GaseousMixture instance
    auto mixture() const -> const GaseousMixture&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
