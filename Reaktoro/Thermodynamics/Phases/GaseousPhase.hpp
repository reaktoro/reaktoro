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
#include <tuple>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivity.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>

namespace Reaktoro {

// Forward declarations
class Phase;

/// Class that defines a gaseous phase
class GaseousPhase : public GaseousMixture
{
public:
    /// Construct a default GaseousPhase instance.
    GaseousPhase();

    /// Construct an GaseousPhase instance with given species.
    explicit GaseousPhase(const std::vector<GaseousSpecies>& species);

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see GaseousActivityFunction
    auto setActivityModel(const std::string& species, const GaseousActivityFunction& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(const std::string& species) -> void;

    /// Set the activity model of CO2(g) to be the one of Duan et al. (2006).
    auto setActivityModelDuanSunCO2() -> void;

    /// Set the activity model of H2O(g) and CO2(g) to be the one of Spycher et al. (2003).
    auto setActivityModelSpycherPruessH2OCO2() -> void;

    /// Set the activity model of H2O(g), CO2(g) and CH4(g) to be the one of Spycher and Reed (1988).
    auto setActivityModelSpycherReedH2OCO2CH4() -> void;

    /// Set the activity model of a gaseous species to be the one of Peng and Robinson (1978).
    /// @param species The name of the gaseous species
    auto setActivityModelPengRobinson(const std::string& species) -> void;

    /// Calculate the concentrations of the gaseous species.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of Pa)
    /// @param n The molar abundance of the species
    /// @return The concentrations of the gaseous species
    auto concentrations(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the activity constants of the gaseous species and their partial derivatives.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of Pa)
    /// @return The activity constants of the gaseous species and their partial derivatives
    auto activityConstants(double T, double P) const -> ThermoVector;

    /// Calculate the activity coefficients of the gaseous species and their partial derivatives.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of Pa)
    /// @param n The molar composition of the gaseous phase
    /// @return The activity coefficients of the gaseous species and their partial derivatives
    auto activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector;

    /// Calculate the activities of the gaseous species and their partial derivatives.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of Pa)
    /// @param n The molar composition of the gaseous phase
    /// @return The activities of the gaseous species and their partial derivatives
    auto activities(double T, double P, const Vector& n) const -> ChemicalVector;

private:
    /// The gaseous activity functions
    std::vector<GaseousActivityFunction> activity_fns;
};

} // namespace Reaktoro
