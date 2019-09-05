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

#include "FluidSpecies.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct FluidSpecies::Impl
{
    // The critical temperature of the fluids (gaseous or liquid) species (in units of K)
    double critical_temperature = 0.0;

    // The critical pressure of the fluids (gaseous or liquid) species (in units of Pa)
    double critical_pressure = 0.0;

    // The acentric factor of the fluids (gaseous or liquid) species
    double acentric_factor = 0.0;

    /// The thermodynamic data of the fluids (gaseous or liquid) species.
    FluidSpeciesThermoData thermo;
};

FluidSpecies::FluidSpecies()
    : pimpl(new Impl())
{}

FluidSpecies::FluidSpecies(const Species& species)
    : Species(species), pimpl(new Impl())
{}

auto FluidSpecies::setCriticalTemperature(double val) -> void
{
    Assert(val > 0.0, "Cannot set the critical temperature of the fluid `" + name() + "`.",
        "The given critical temperature `" + std::to_string(val) + "` is not positive.");
    pimpl->critical_temperature = val;
}

auto FluidSpecies::setCriticalPressure(double val) -> void
{
    Assert(val > 0.0, "Cannot set the critical pressure of the fluid `" + name() + "`.",
        "The given critical pressure `" + std::to_string(val) + "` is not positive.");
    pimpl->critical_pressure = val;
}

auto FluidSpecies::setAcentricFactor(double val) -> void
{
    pimpl->acentric_factor = val;
}

auto FluidSpecies::setThermoData(const FluidSpeciesThermoData& thermo) -> void
{
    pimpl->thermo = thermo;
}

auto FluidSpecies::criticalTemperature() const -> double
{
    return pimpl->critical_temperature;
}

auto FluidSpecies::criticalPressure() const -> double
{
    return pimpl->critical_pressure;
}

auto FluidSpecies::acentricFactor() const -> double
{
    return pimpl->acentric_factor;
}

auto FluidSpecies::thermoData() const -> const FluidSpeciesThermoData&
{
    return pimpl->thermo;
}

} // namespace Reaktoro
