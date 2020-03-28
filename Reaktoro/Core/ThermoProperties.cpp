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

#include "ThermoProperties.hpp"

namespace Reaktoro {

ThermoProperties::ThermoProperties()
{}

ThermoProperties::ThermoProperties(const ChemicalSystem& system)
: system(system), num_species(system.numSpecies()), num_phases(system.numPhases()),
  tres(num_phases, num_species), T(298.15), P(1e-5)
{}

auto ThermoProperties::update(double T_, double P_) -> void
{
    // Update both temperature and pressure
    T = T_;
    P = P_;

    // Update the thermodynamic properties of each phase
    Index ispecies = 0;
    for(Index iphase = 0; iphase < num_phases; ++iphase)
    {
        const auto nspecies = system.numSpeciesInPhase(iphase);
        auto tp = tres.phaseProperties(iphase, ispecies, nspecies);
        system.phase(iphase).thermoModel()(tp, T, P);
        ispecies += nspecies;
    }
}

auto ThermoProperties::temperature() const -> Temperature
{
    return T;
}

auto ThermoProperties::pressure() const -> Pressure
{
    return P;
}

auto ThermoProperties::standardPartialMolarGibbsEnergies() const -> VectorXr
{
    return tres.standardPartialMolarGibbsEnergies();
}

auto ThermoProperties::standardPartialMolarEnthalpies() const -> VectorXr
{
    return tres.standardPartialMolarEnthalpies();
}

auto ThermoProperties::standardPartialMolarVolumes() const -> VectorXr
{
    return tres.standardPartialMolarVolumes();
}

auto ThermoProperties::standardPartialMolarEntropies() const -> VectorXr
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& H = standardPartialMolarEnthalpies();
    return (H - G)/T;
}

auto ThermoProperties::standardPartialMolarInternalEnergies() const -> VectorXr
{
    const auto& H = standardPartialMolarEnthalpies();
    const auto& V = standardPartialMolarVolumes();
    return H - P*V;
}

auto ThermoProperties::standardPartialMolarHelmholtzEnergies() const -> VectorXr
{
    const auto& G = standardPartialMolarGibbsEnergies();
    const auto& V = standardPartialMolarVolumes();
    return G - P*V;
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstP() const -> VectorXr
{
    return tres.standardPartialMolarHeatCapacitiesConstP();
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstV() const -> VectorXr
{
    return tres.standardPartialMolarHeatCapacitiesConstV();
}

} // namespace Reaktoro
