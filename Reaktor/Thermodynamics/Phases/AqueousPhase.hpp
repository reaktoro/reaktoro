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

// C++ includes
#include <map>
#include <string>
#include <tuple>
#include <vector>

// Reaktor includes
#include <Reaktor/Activity/AqueousActivity.hpp>
#include <Reaktor/Mixtures/AqueousMixture.hpp>

namespace Reaktor {

// Reaktor forward declarations
class Phase;

/// Class that defines an aqueous phase
class AqueousPhase : public AqueousMixture
{
public:
    /// Construct a default AqueousPhase instance.
    AqueousPhase();

    /// Construct an AqueousPhase instance with given species.
    explicit AqueousPhase(const std::vector<AqueousSpecies>& species);

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see AqueousActivity
    auto setActivityModel(const std::string& species, const AqueousActivity& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(const std::string& species) -> void;

    /// Set the activity model of the species to be the Setschenow one.
    /// @param species The name of species to have its activity model set
    /// @param b The Setschenow constant
    auto setActivityModelSetschenow(const std::string& species, double b) -> void;

    /// Set the activity model of CO2(aq) to be the one of Duan and Sun (2003).
    auto setActivityModelDuanSunCO2() -> void;

    /// Set the activity model of CO2(aq) to be the one of Drummond (1981).
    auto setActivityModelDrummondCO2() -> void;

    /// Set the activity model of CO2(aq) to be the one of Rumpf et al. (1994).
    auto setActivityModelRumpfCO2() -> void;

    /// Set the activity model of H2O(l) to be the one of Helgeson et al. (1981).
    auto setActivityModelHKFWater() -> void;

    /// Set the activity models of charged species to be the one of Helgeson et al. (1981).
    auto setActivityModelHKFChargedSpecies() -> void;

    /// Set the activity model of H2O(l) to the Pitzer's model of Harvie et al. (1984).
    auto setActivityModelPitzerWater() -> void;

    /// Set the activity model of the charged species to the Pitzer's model of Harvie et al. (1984).
    auto setActivityModelPitzerChargedSpecies() -> void;

    /// Set the activity model of a neutral species to the Pitzer's model of Harvie et al. (1984).
    /// @param species The name of the neutral species
    auto setActivityModelPitzerNeutralSpecies(const std::string& species) -> void;

    /// Calculate the parameters for the aqueous activity calculation.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of bar)
    /// @param n The molar composition of the aqueous phase
    /// @return The parameters to be used for the aqueous activity calculation
    auto params(double T, double P, const Vector& n) const -> AqueousActivityParams;

    /// Calculate the concentrations of the aqueous species.
    /// @param n The molar abundance of the species
    /// @return The concentrations of the aqueous species
    auto concentrations(const Vector& n) const -> Vector;

    /// Calculate the activities of the aqueous species and its molar derivatives.
    /// @param T The temperature used for the calculation (in units of K)
    /// @param P The pressure used for the calculation (in units of bar)
    /// @param n The molar composition of the aqueous phase
    /// @return The activities of the aqueous species and their molar derivatives
    auto activities(double T, double P, const Vector& n) const -> PartialVector;

private:
    /// The aqueous activity functions
    std::vector<AqueousActivity> activities$;
};

/// Creates a Phase instance from an AqueousPhase instance
/// @param phase The AqueousPhase instance
/// @return A Phase instance created from the given aqueous phase
auto createPhase(const AqueousPhase& phase) -> Phase;

} // namespace Reaktor