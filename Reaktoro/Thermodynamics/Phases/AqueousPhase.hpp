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
#include <Reaktoro/Thermodynamics/Activity/AqueousActivity.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;

/// Class that defines an aqueous phase
class AqueousPhase : public Phase
{
public:
    /// Construct a default AqueousPhase instance.
    AqueousPhase();

    /// Contruct a copy of an AqueousPhase instance
    AqueousPhase(const AqueousPhase& other);

    /// Construct an AqueousPhase instance with given aqueous mixture.
    explicit AqueousPhase(const AqueousMixture& mixture);

    /// Destroy this instance
    virtual ~AqueousPhase();

    /// Assign an AqueousPhase instance to this
    auto operator=(AqueousPhase other) -> AqueousPhase&;

    /// Set the chemical model of the phase with the HKF equation of state.
    auto setChemicalModelHKF() -> void;

    /// Set the chemical model of the phase with the Pitzer equation of state.
    auto setChemicalModelPitzer() -> void;

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see AqueousActivityFunction
    auto setActivityModel(std::string species, const AqueousActivityFunction& activity) -> void;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(std::string species) -> void;

    /// Set the activity model of the species to be the Setschenow one.
    /// @param species The name of species to have its activity model set
    /// @param b The Setschenow constant
    auto setActivityModelSetschenow(std::string species, double b) -> void;

    /// Set the activity model of CO2(aq) to be the one of Duan and Sun (2003).
    auto setActivityModelDuanSunCO2() -> void;

    /// Set the activity model of CO2(aq) to be the one of Drummond (1981).
    auto setActivityModelDrummondCO2() -> void;

    /// Set the activity model of CO2(aq) to be the one of Rumpf et al. (1994).
    auto setActivityModelRumpfCO2() -> void;

    /// Return the AqueousMixture instance
    auto mixture() const -> const AqueousMixture&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
