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
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;
class DebyeHuckelParams;

/// A type used to describe an aqueous phase.
class AqueousPhase : public Phase
{
public:
    /// Construct a default AqueousPhase instance.
    AqueousPhase();

    /// Construct an AqueousPhase instance with given aqueous mixture.
    explicit AqueousPhase(const AqueousMixture& mixture);

    /// Set the temperature and pressure interpolation points for calculation of water density and water dielectric constant.
    /// Use this method if temperature-pressure interpolation should be used for the calculation of water density and
    /// water dielectric constant. This should be done if the cost of the analytical calculation of these properties
    /// is prohibitive for your application.
    /// @param temperatures The temperature points (in units of K)
    /// @param pressures The pressure points (in units of Pa)
    auto setInterpolationPoints(const std::vector<double>& temperatures, const std::vector<double>& pressures) -> AqueousPhase&;

    /// Set the chemical model of the phase with the ideal aqueous solution equation of state.
    auto setChemicalModelIdeal() -> AqueousPhase&;

    /// Set the chemical model of the phase with the Debye-Huckel equation of state.
    /// {
    auto setChemicalModelDebyeHuckel() -> AqueousPhase&;
    auto setChemicalModelDebyeHuckel(const DebyeHuckelParams& params) -> AqueousPhase&;
    /// }

    /// Set the chemical model of the phase with the HKF equation of state.
    auto setChemicalModelHKF() -> AqueousPhase&;

    /// Set the chemical model of the phase with the Pitzer equation of state.
    /// Uses the Pitzer equation of state described in:
    /// *Harvie, C.E., Møller, N., Weare, J.H. (1984). The prediction of mineral
    /// solubilities in natural waters: The Na-K-Mg-Ca-H-Cl-SO4-OH-HCO3-CO3-CO2-H2O
    /// system to high ionic strengths at 25°C. Geochimica et Cosmochimica Acta, 48(4), 723–751*.
    auto setChemicalModelPitzerHMW() -> AqueousPhase&;

    /// Set the activity model of a species.
    /// @param species The name of the species
    /// @param activity The activity function
    /// @see AqueousActivityModel
    auto setActivityModel(std::string species, const AqueousActivityModel& activity) -> AqueousPhase&;

    /// Set the activity model of the species to be the ideal one.
    /// @param species The name of species to have its activity model set
    auto setActivityModelIdeal(std::string species) -> AqueousPhase&;

    /// Set the activity model of the species to be the Setschenow one.
    /// @param species The name of species to have its activity model set
    /// @param b The Setschenow constant
    auto setActivityModelSetschenow(std::string species, double b) -> AqueousPhase&;

    /// Set the activity model of CO2(aq) to be the one of Duan and Sun (2003).
    auto setActivityModelDuanSunCO2() -> AqueousPhase&;

    /// Set the activity model of CO2(aq) to be the one of Drummond (1981).
    auto setActivityModelDrummondCO2() -> AqueousPhase&;

    /// Set the activity model of CO2(aq) to be the one of Rumpf et al. (1994).
    auto setActivityModelRumpfCO2() -> AqueousPhase&;

    /// Return the AqueousMixture instance
    auto mixture() const -> const AqueousMixture&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
